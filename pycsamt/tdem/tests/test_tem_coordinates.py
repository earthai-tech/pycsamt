"""Tests for TEM profile/point coordinate readers."""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.tdem import (
    TEMCoordinateTable,
    read_tem_coordinates,
    read_temavg_survey,
)
from pycsamt.tdem.io import (
    read_tem_coordinates as read_tem_coordinates_io,
)

DATA_DIR = Path(__file__).parents[3] / "data" / "TEMAVG" / "JIANGSU"
AVG_FILE = DATA_DIR / "TEM100.AVG"


pytestmark = pytest.mark.skipif(
    not AVG_FILE.exists(),
    reason="TEMAVG sample data not available",
)


def _write_coordinate_csv(path: Path) -> None:
    """Create a small coordinate table in the field layout."""
    path.write_text(
        "\n".join(
            [
                "Profile,point,Gauss Coordinate,,Relative coordinate,,H(m),",
                ",,X(m),Y(m),X(m),Y(m),,",
                "100,100,4291789.7679,19510112.9006," "100.0034,100.0151,1102.9537,",
                "100,120,4291809.1,19510120.2,120.0,100.5,1103.1,road",
            ]
        ),
        encoding="utf-8",
    )


def test_read_tem_coordinates_csv(tmp_path):
    """Coordinate CSV files should be parsed by profile and point."""
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)

    coords = read_tem_coordinates(coord_path)
    first = coords.get(100, 100)
    second = coords.get(100, 120)

    assert isinstance(coords, TEMCoordinateTable)
    assert coords.n_points == 2
    assert coords.profiles == [100.0]
    assert coords.points == [100.0, 120.0]
    assert first is not None
    assert first.gauss_x == pytest.approx(4291789.7679)
    assert first.gauss_y == pytest.approx(19510112.9006)
    assert first.x == pytest.approx(100.0034)
    assert first.y == pytest.approx(100.0151)
    assert first.elevation == pytest.approx(1102.9537)
    assert second is not None
    assert second.remark == "road"


def test_read_tem_coordinates_io_wrapper(tmp_path):
    """The public IO module should expose the coordinate reader."""
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)

    coords = read_tem_coordinates_io(coord_path)

    assert coords.n_points == 2
    assert coords.get(100, 100) is not None


def test_survey_records_are_enriched_from_coordinate_file(tmp_path):
    """Explicit coordinates should be attached to matching AVG rows."""
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)

    survey = read_temavg_survey(DATA_DIR, coordinate_file=coord_path)
    coord = survey.coordinate_for(100, 100)
    first = survey.to_records()[0]

    assert survey.coordinates is not None
    assert coord is not None
    assert first["profile"] == pytest.approx(100.0)
    assert first["station"] == pytest.approx(100.0)
    assert first["coord_profile"] == pytest.approx(100.0)
    assert first["coord_point"] == pytest.approx(100.0)
    assert first["x"] == pytest.approx(100.0034)
    assert first["y"] == pytest.approx(100.0151)
    assert first["elevation"] == pytest.approx(1102.9537)


def test_survey_soundings_are_enriched_from_coordinate_file(tmp_path):
    """Generated soundings should carry matching station coordinates."""
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)
    survey = read_temavg_survey(DATA_DIR, coordinate_file=coord_path)

    soundings = survey.to_soundings(stems=["TEM100"])
    first = soundings[0]

    assert first.station_name == "TEM100_100"
    assert first.x == pytest.approx(100.0034)
    assert first.y == pytest.approx(100.0151)
    assert first.elevation == pytest.approx(1102.9537)
