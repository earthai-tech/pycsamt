# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared PCPT (points / targets) app-adapter tests."""

from __future__ import annotations

import base64
import json

import pytest

from pycsamt.app._points import (
    _hover,
    decode_points_upload,
    point_map_markers,
    point_scene_traces,
    pointset_from_store,
)
from pycsamt.format import Point, PointSet, point_set_to_dict
from pycsamt.map.geometry import survey_frame


def _set() -> PointSet:
    return PointSet(
        document_id="pcpt:t",
        created_at="2026-09-03T00:00:00Z",
        created_by="pytest",
        crs="EPSG:4326",
        points=[
            Point("T1", "Target 1", longitude=103.001, latitude=22.001,
                  depth_top=100.0, depth_bottom=350.0, kind="target"),
            Point("T2", "Target 2", longitude=103.002, latitude=22.002),
        ],
    )


def _upload(point_set: PointSet) -> str:
    raw = json.dumps(point_set_to_dict(point_set)).encode()
    return "data:application/json;base64," + base64.b64encode(raw).decode()


def test_decode_upload_round_trips():
    store = decode_points_upload(_upload(_set()), "targets.pcpt.json")
    assert store["n_points"] == 2
    restored = pointset_from_store(store)
    assert [p.id for p in restored.points] == ["T1", "T2"]


def test_map_markers_one_trace_coloured_by_kind():
    traces = point_map_markers(_set())
    assert len(traces) == 1
    colours = list(traces[0].marker.color)
    assert colours[0] == "#ef4444"  # target


def test_scene_traces_place_points_and_draw_depth_stems():
    frame = survey_frame(
        [22.0, 22.001, 22.002],
        [103.0, 103.001, 103.002],
        ["L", "L", "L"],
    )
    traces = point_scene_traces(_set(), frame, datum="zero")
    kinds = [t.type for t in traces]
    assert "scatter3d" in kinds
    # T1 has a depth window -> one stem line in addition to the marker trace
    assert len(traces) == 2


def test_scene_traces_empty_without_frame():
    assert point_scene_traces(_set(), None) == []


# ── decode_points_upload: base64 / branch errors ──────────────────────────


def test_decode_upload_rejects_missing_comma():
    with pytest.raises(ValueError, match="invalid PCPT upload payload"):
        decode_points_upload("not-a-data-url", "targets.pcpt.json")


def test_decode_upload_rejects_bad_base64():
    with pytest.raises(ValueError, match="invalid base64 PCPT upload"):
        decode_points_upload(
            "data:application/json;base64,not*valid*base64!!", "t.pcpt.json"
        )


def test_decode_upload_rejects_oversized_payload(monkeypatch):
    import pycsamt.app._points as points_mod

    monkeypatch.setattr(points_mod, "_MAX_BYTES", 4)
    raw = base64.b64encode(b"way too many bytes").decode()
    with pytest.raises(ValueError, match="exceeds"):
        decode_points_upload(
            "data:application/json;base64," + raw, "t.pcpt.json"
        )


def test_decode_upload_from_csv(tmp_path):
    csv_text = (
        "id,name,longitude,latitude\n"
        "T1,Target 1,103.001,22.001\n"
        "T2,Target 2,103.002,22.002\n"
    )
    payload = base64.b64encode(csv_text.encode()).decode()
    store = decode_points_upload(
        "data:text/csv;base64," + payload, "targets.csv"
    )
    assert store["n_points"] == 2
    restored = pointset_from_store(store)
    assert [p.id for p in restored.points] == ["T1", "T2"]


def test_decode_upload_from_xlsx(monkeypatch):
    import pycsamt.app._points as points_mod

    def _fake_points_from_xlsx(raw):
        return _set()

    monkeypatch.setattr(
        "pycsamt.format.points_from_xlsx", _fake_points_from_xlsx
    )
    payload = base64.b64encode(b"fake-xlsx-bytes").decode()
    store = decode_points_upload(
        "data:application/octet-stream;base64," + payload, "targets.xlsx"
    )
    assert store["n_points"] == 2


# ── pointset_from_store edge cases ─────────────────────────────────────────


def test_pointset_from_store_none_returns_none():
    assert pointset_from_store(None) is None


def test_pointset_from_store_empty_dict_returns_none():
    assert pointset_from_store({}) is None


def test_pointset_from_store_rejects_non_dict_payload():
    with pytest.raises(ValueError, match="invalid PCPT application store"):
        pointset_from_store({"document": "not-a-dict"})


# ── point_map_markers edge cases ───────────────────────────────────────────


def test_map_markers_empty_without_store():
    assert point_map_markers(None) == []


def test_map_markers_empty_when_no_locatable_points():
    point_set = PointSet(
        document_id="pcpt:noloc",
        created_at="2026-09-03T00:00:00Z",
        created_by="pytest",
        crs="",
        points=[Point("P1", "No location", x=100.0, y=200.0)],
    )
    assert point_map_markers(point_set) == []


def test_map_markers_accepts_store_dict():
    store = {"document": point_set_to_dict(_set())}
    traces = point_map_markers(store)
    assert len(traces) == 1


# ── point_scene_traces datum handling ──────────────────────────────────────


def test_scene_traces_surface_datum_uses_surface_callable():
    frame = survey_frame(
        [22.0, 22.001, 22.002],
        [103.0, 103.001, 103.002],
        ["L", "L", "L"],
    )
    traces = point_scene_traces(
        _set(), frame, surface=lambda x: 50.0, datum="surface"
    )
    marker = traces[0]
    # T1 has depth_top=100 -> top = 50 - 100; T2 has no depth -> top = 50
    assert sorted(marker.z) == [pytest.approx(-50.0), pytest.approx(50.0)]


def test_scene_traces_collar_z_datum_uses_point_z():
    point_set = PointSet(
        document_id="pcpt:z",
        created_at="2026-09-03T00:00:00Z",
        created_by="pytest",
        crs="EPSG:4326",
        points=[
            Point(
                "T1", "Target 1", longitude=103.001, latitude=22.001,
                z=250.0, kind="target",
            ),
        ],
    )
    frame = survey_frame(
        [22.0, 22.001, 22.002],
        [103.0, 103.001, 103.002],
        ["L", "L", "L"],
    )
    traces = point_scene_traces(point_set, frame, datum="collar_z")
    assert traces[0].z[0] == pytest.approx(250.0)


def test_scene_traces_empty_without_points():
    frame = survey_frame(
        [22.0, 22.001, 22.002],
        [103.0, 103.001, 103.002],
        ["L", "L", "L"],
    )
    point_set = PointSet(
        document_id="pcpt:empty",
        created_at="2026-09-03T00:00:00Z",
        created_by="pytest",
        crs="",
        points=[],
    )
    assert point_scene_traces(point_set, frame) == []


def test_scene_traces_accepts_store_dict():
    frame = survey_frame(
        [22.0, 22.001, 22.002],
        [103.0, 103.001, 103.002],
        ["L", "L", "L"],
    )
    store = {"document": point_set_to_dict(_set())}
    traces = point_scene_traces(store, frame, datum="zero")
    assert len(traces) == 2


# ── _hover ──────────────────────────────────────────────────────────────────


def test_hover_none_point_returns_pid():
    assert _hover(None, "P1") == "P1"


def test_hover_includes_note_and_single_depth():
    point = Point(
        "T1", "Target 1", longitude=1.0, latitude=1.0,
        depth_top=10.0, note="drill here",
    )
    text = _hover(point, "T1")
    assert "depth 10 m" in text
    assert "drill here" in text
