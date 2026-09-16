# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for the ``pycsamt map`` command group."""

from __future__ import annotations

import csv
import json
import math
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from pycsamt.cli import main
from pycsamt.cli.commands.map._base import (
    _finite_float,
    _format_number,
    _format_rows,
    _is_empty,
    _site_summary,
    _station_rows,
    _station_rows_from_map_api,
)
from pycsamt.cli.tests.conftest import make_fake_sites


def test_map_group_help(runner: CliRunner) -> None:
    result = runner.invoke(main, ["map", "--help"])

    assert result.exit_code == 0
    assert "stations" in result.output
    assert "plot" in result.output


def test_map_registered_on_root_help(runner: CliRunner) -> None:
    result = runner.invoke(main, ["--help"])

    assert result.exit_code == 0
    assert "map" in result.output


def test_map_stations_help(runner: CliRunner) -> None:
    result = runner.invoke(main, ["map", "stations", "--help"])

    assert result.exit_code == 0
    assert "--drop-missing" in result.output
    assert "--format" in result.output
    assert "--output" in result.output


def test_map_plot_help(runner: CliRunner) -> None:
    result = runner.invoke(main, ["map", "plot", "--help"])

    assert result.exit_code == 0
    assert "--output" in result.output
    assert "--no-label" in result.output
    assert "--dpi" in result.output


def test_map_stations_json(runner: CliRunner, monkeypatch) -> None:
    fake = make_fake_sites(3)
    monkeypatch.setattr(
        "pycsamt.cli.commands.map.stations._get_sites",
        lambda *a, **kw: fake,
    )

    result = runner.invoke(
        main,
        ["map", "stations", "--survey", ".", "--format", "json"],
    )

    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data) == 3
    assert data[0]["station"] == "S01"
    assert data[0]["lat"] == 0.0
    assert data[0]["lon"] == 100.0


def test_map_stations_csv(runner: CliRunner, monkeypatch) -> None:
    fake = make_fake_sites(2)
    monkeypatch.setattr(
        "pycsamt.cli.commands.map.stations._get_sites",
        lambda *a, **kw: fake,
    )

    result = runner.invoke(
        main,
        ["map", "stations", "--survey", ".", "--format", "csv"],
    )

    assert result.exit_code == 0
    rows = list(csv.DictReader(StringIO(result.output)))
    assert len(rows) == 2
    assert rows[1]["station"] == "S02"
    assert rows[1]["lat"] == "1.0"
    assert rows[1]["lon"] == "101.0"


def test_map_stations_text_default_format(runner: CliRunner, monkeypatch) -> None:
    fake = make_fake_sites(2)
    monkeypatch.setattr(
        "pycsamt.cli.commands.map.stations._get_sites",
        lambda *a, **kw: fake,
    )

    result = runner.invoke(main, ["map", "stations", "--survey", "."])

    assert result.exit_code == 0
    assert "station" in result.output
    assert "S01" in result.output
    assert "S02" in result.output


def test_map_stations_top_limits_rows(runner: CliRunner, monkeypatch) -> None:
    fake = make_fake_sites(5)
    monkeypatch.setattr(
        "pycsamt.cli.commands.map.stations._get_sites",
        lambda *a, **kw: fake,
    )

    result = runner.invoke(
        main,
        ["map", "stations", "--survey", ".", "--format", "json", "--top", "2"],
    )

    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data) == 2


def test_map_stations_output_file(
    runner: CliRunner, monkeypatch, tmp_path: Path
) -> None:
    fake = make_fake_sites(2)
    monkeypatch.setattr(
        "pycsamt.cli.commands.map.stations._get_sites",
        lambda *a, **kw: fake,
    )
    out = tmp_path / "coords.csv"

    result = runner.invoke(
        main,
        [
            "map",
            "stations",
            "--survey",
            ".",
            "--format",
            "csv",
            "--output",
            str(out),
        ],
    )

    assert result.exit_code == 0
    assert "Wrote 2 station row(s)" in result.output
    assert out.exists()
    assert "S01" in out.read_text(encoding="utf-8")


def test_map_stations_output_write_error(
    runner: CliRunner, monkeypatch, tmp_path: Path
) -> None:
    fake = make_fake_sites(1)
    monkeypatch.setattr(
        "pycsamt.cli.commands.map.stations._get_sites",
        lambda *a, **kw: fake,
    )
    # A path whose parent cannot be created (colliding with an existing file)
    bogus_parent = tmp_path / "not_a_dir"
    bogus_parent.write_text("x")
    out = bogus_parent / "coords.csv"

    result = runner.invoke(
        main,
        ["map", "stations", "--survey", ".", "--output", str(out)],
    )

    assert result.exit_code == 1
    assert "Error writing" in result.output


def test_map_stations_real_survey(runner: CliRunner, site_edi_dir: Path) -> None:
    """Exercise ``_get_sites`` -> ``resolve_survey`` with real EDI data."""
    result = runner.invoke(
        main, ["map", "stations", str(site_edi_dir), "--format", "json"]
    )
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data) >= 1


# ---------------------------------------------------------------------------
# pycsamt map plot
# ---------------------------------------------------------------------------


def test_map_plot_real_survey_saves_png(
    runner: CliRunner, site_edi_dir: Path, tmp_path: Path
) -> None:
    out = tmp_path / "map.png"
    result = runner.invoke(
        main,
        [
            "map",
            "plot",
            str(site_edi_dir),
            "--output",
            str(out),
            "--no-label",
            "--dpi",
            "75",
            "--title",
            "MyMap",
        ],
    )
    assert result.exit_code == 0, result.output
    assert out.exists()
    assert "Saved station map" in result.output


def test_map_plot_show_flag(
    runner: CliRunner, site_edi_dir: Path, tmp_path: Path
) -> None:
    out = tmp_path / "map.png"
    with patch("matplotlib.pyplot.show") as mock_show:
        result = runner.invoke(
            main,
            ["map", "plot", str(site_edi_dir), "--output", str(out), "--show"],
        )
    assert result.exit_code == 0, result.output
    mock_show.assert_called_once()


def test_map_plot_no_stations_found(runner: CliRunner) -> None:
    class _EmptySites:
        def __len__(self) -> int:
            return 0

        def __iter__(self):
            return iter([])

    with patch(
        "pycsamt.cli.commands.map.plot._get_sites",
        return_value=_EmptySites(),
    ):
        result = runner.invoke(main, ["map", "plot", "--survey", "."])

    assert result.exit_code != 0
    assert "No stations with finite latitude and longitude" in result.output


def test_map_plot_import_error(runner: CliRunner) -> None:
    fake = make_fake_sites(2)
    with (
        patch("pycsamt.cli.commands.map.plot._get_sites", return_value=fake),
        patch(
            "pycsamt.map.ensure_map_data",
            side_effect=ImportError("optional dependency missing"),
        ),
    ):
        result = runner.invoke(main, ["map", "plot", "--survey", "."])

    assert result.exit_code == 1
    assert "optional dependency missing" in result.output


def test_map_plot_generic_error(runner: CliRunner) -> None:
    fake = make_fake_sites(2)
    with (
        patch("pycsamt.cli.commands.map.plot._get_sites", return_value=fake),
        patch("pycsamt.map.ensure_map_data", side_effect=RuntimeError("boom")),
    ):
        result = runner.invoke(main, ["map", "plot", "--survey", "."])

    assert result.exit_code == 1
    assert "Error building station map: boom" in result.output


# ---------------------------------------------------------------------------
# _base helper unit tests
# ---------------------------------------------------------------------------


def test_finite_float_valid() -> None:
    assert _finite_float("3.5") == 3.5


def test_finite_float_non_numeric_returns_none() -> None:
    assert _finite_float("abc") is None
    assert _finite_float(None) is None


def test_finite_float_nan_inf_returns_none() -> None:
    assert _finite_float(math.nan) is None
    assert _finite_float(math.inf) is None


def test_format_number_none_returns_empty() -> None:
    assert _format_number(None, 2) == ""


def test_format_number_non_numeric_falls_back_to_str() -> None:
    assert _format_number("n/a", 2) == "n/a"


def test_format_number_formats_float() -> None:
    assert _format_number(3.14159, 2) == "3.14"


def test_format_rows_empty_text() -> None:
    assert _format_rows([], "text") == "No stations found."


def test_is_empty_non_sized_object_returns_false() -> None:
    assert _is_empty(object()) is False


def test_site_summary_summary_method_used() -> None:
    class _S:
        def summary(self):
            return {"name": "X1", "lat": 1.0, "lon": 2.0, "elev": 3.0, "nfreq": 5}

    out = _site_summary(_S())
    assert out["name"] == "X1"


def test_site_summary_coords_exception_falls_back() -> None:
    class _Broken:
        @property
        def coords(self):
            raise RuntimeError("no coords")

    out = _site_summary(_Broken())
    assert out["lat"] is None
    assert out["lon"] is None


def test_station_rows_from_map_api_exception_returns_none() -> None:
    fake = make_fake_sites(2)
    with patch("pycsamt.map.ensure_map_data", side_effect=RuntimeError("x")):
        assert _station_rows_from_map_api(fake, drop_missing=False) is None


def test_station_rows_drop_missing_manual_path() -> None:
    fake = make_fake_sites(3)
    rows = _station_rows(fake, drop_missing=True)
    assert len(rows) == 3
    for row in rows:
        assert row["lat"] is not None and row["lon"] is not None
