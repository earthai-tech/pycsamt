"""Tests for the PCPT points / targets mini-format."""

from __future__ import annotations

import io

import pytest

from pycsamt.format import (
    Point,
    PointSet,
    PointSetValidationError,
    point_set_from_dict,
    point_set_to_dict,
    points_from_csv,
    points_from_xlsx,
    read_points,
    write_points,
)


def _set() -> PointSet:
    return PointSet(
        document_id="pcpt:test",
        created_at="2026-09-03T00:00:00Z",
        created_by="pytest",
        crs="EPSG:32648",
        points=[
            Point("T1", "Target 1", x=500.0, y=1000.0, z=300.0,
                  depth_top=150.0, depth_bottom=400.0, kind="target"),
            Point("T2", "Target 2", x=560.0, y=1040.0),
        ],
    )


def test_json_round_trip(tmp_path):
    original = _set()
    path = tmp_path / "targets.pcpt.json"
    write_points(original, path)
    restored = read_points(path)
    assert [p.id for p in restored.points] == ["T1", "T2"]
    assert restored.points[0].depth_bottom == 400.0
    assert restored.crs == "EPSG:32648"


def test_to_dict_omits_empty_optionals_but_keeps_id():
    payload = point_set_to_dict(_set())
    assert payload["points"][1] == {"id": "T2", "name": "Target 2",
                                    "x": 560.0, "y": 1040.0,
                                    "kind": "poi", "symbol": "circle"}


def test_validation_rejects_point_without_location():
    bad = PointSet(
        document_id="pcpt:bad",
        created_at="2026-09-03T00:00:00Z",
        created_by="pytest",
        crs="EPSG:4326",
        points=[Point("X", "no location")],
    )
    with pytest.raises(PointSetValidationError):
        bad.validate()


def test_lonlat_passthrough_and_reprojection():
    wgs = PointSet(
        document_id="pcpt:w",
        created_at="2026-09-03T00:00:00Z",
        created_by="pytest",
        crs="EPSG:4326",
        points=[Point("A", longitude=103.1, latitude=22.2)],
    )
    assert wgs.lonlat()["A"] == (22.2, 103.1)


def test_points_from_csv_with_aliases(tmp_path):
    csv_path = tmp_path / "targets.csv"
    csv_path.write_text(
        "Target,Longitude,Latitude,Depth from,Depth to,Note\n"
        "ZK-P1,103.10,22.20,200,450,porphyry\n"
        "ZK-P2,103.11,22.21,,,skarn\n",
        encoding="utf-8",
    )
    point_set = points_from_csv(csv_path)
    assert point_set.crs == "EPSG:4326"
    assert [p.id for p in point_set.points] == ["ZK-P1", "ZK-P2"]
    assert point_set.points[0].depth_bottom == 450.0
    assert point_set.points[0].note == "porphyry"


def test_points_from_xlsx_projected(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["Target", "Easting", "Northing", "Elev"])
    ws.append(["P1", 350000, 3000000, 1200])
    ws.append(["P2", 350120, 3000050, 1210])
    buffer = io.BytesIO()
    wb.save(buffer)
    point_set = points_from_xlsx(buffer.getvalue(), crs="EPSG:32648")
    assert point_set.crs == "EPSG:32648"
    assert point_set.points[1].x == 350120.0


def test_from_dict_rejects_non_object():
    with pytest.raises(PointSetValidationError):
        point_set_from_dict([1, 2, 3])
