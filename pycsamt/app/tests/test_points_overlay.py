# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared PCPT (points / targets) app-adapter tests."""

from __future__ import annotations

import base64
import json

from pycsamt.app._points import (
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
