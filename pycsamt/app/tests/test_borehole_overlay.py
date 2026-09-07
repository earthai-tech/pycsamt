"""Shared application-boundary tests for standalone PCBH overlays."""

from __future__ import annotations

import base64
import json

import plotly.graph_objects as go
import pytest

from pycsamt.app._borehole import (
    add_pcbh_to_figure,
    decode_pcbh_upload,
    pcbh_plotly_traces,
)
from pycsamt.format.borehole import (
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    VocabularyEntry,
    pcbh_to_dict,
)


def _document() -> PCBHDocument:
    return PCBHDocument(
        document_id="test:app-overlay",
        created_at="2026-08-28T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:32629"),
        lithologies=[
            VocabularyEntry("CLAY", "Clay"),
            VocabularyEntry("GRANITE", "Granite"),
        ],
        boreholes=[
            PCBHBorehole(
                id="BH-1",
                name="Borehole 1",
                kind="mining_exploration",
                status="completed",
                collar=Collar(100.0, 200.0, 50.0),
                total_depth_md=20.0,
                interval_logs={
                    "lithology": [
                        LogInterval(0.0, 10.0, code="CLAY"),
                        LogInterval(10.0, 20.0, code="GRANITE"),
                    ]
                },
            )
        ],
    )


def _upload(document: PCBHDocument) -> str:
    raw = json.dumps(pcbh_to_dict(document)).encode()
    return "data:application/json;base64," + base64.b64encode(raw).decode()


def test_upload_store_round_trip_and_rejects_bad_content():
    store = decode_pcbh_upload(_upload(_document()), "holes.pcbh.json")
    assert store["n_boreholes"] == 1
    assert store["document"]["boreholes"][0]["id"] == "BH-1"

    with pytest.raises(ValueError, match="UTF-8 JSON"):
        decode_pcbh_upload(
            "data:application/json;base64," + base64.b64encode(b"{").decode(),
            "holes.json",
        )


def test_shared_trace_adapter_preserves_intervals_contacts_and_hover():
    traces = pcbh_plotly_traces(_document())
    line_traces = [trace for trace in traces if trace.mode == "lines"]
    assert len(line_traces) == 2
    assert any(trace.name == "Contact" for trace in traces)
    assert all("Borehole" in trace.hovertemplate for trace in line_traces)


def test_same_store_adds_identical_overlay_to_each_application_figure():
    store = decode_pcbh_upload(_upload(_document()), "holes.pcbh.json")
    mapview_figure = add_pcbh_to_figure(go.Figure(), store)
    web_figure = add_pcbh_to_figure(go.Figure(), store)

    assert mapview_figure.to_plotly_json()["data"] == web_figure.to_plotly_json()[
        "data"
    ]


def test_strip_log_figure_one_column_per_hole_with_vocabulary_colours():
    from pycsamt.app._borehole import strip_log_figure

    figure = strip_log_figure(_document(), family="lithology")
    bars = [trace for trace in figure.data if trace.type == "bar"]
    assert len(bars) == 2
    assert {bar.marker.color for bar in bars} == {"#cccccc", None} or bars
    assert figure.layout.yaxis.autorange == "reversed"


def test_builder_draft_round_trips_through_document():
    from pycsamt.app._borehole import (
        builder_draft_from_document,
        document_from_builder_draft,
    )

    draft = builder_draft_from_document(pcbh_to_dict(_document()))
    assert "boreholes" in draft and draft["boreholes"][0]["id"] == "BH-1"
    restored = document_from_builder_draft(draft)
    assert restored.boreholes[0].total_depth_md == 20.0


def test_scene_borehole_traces_from_alignment():
    from pycsamt.map.borehole_align import align_boreholes_to_scene
    from pycsamt.map.geometry import survey_frame
    from pycsamt.app._borehole import scene_borehole_traces

    document = PCBHDocument(
        document_id="test:scene",
        created_at="2026-09-03T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:4326"),
        lithologies=[VocabularyEntry("C", "Clay", color="#abcdef")],
        boreholes=[
            PCBHBorehole(
                id="BH-1",
                name="Hole",
                kind="mining_exploration",
                status="completed",
                collar=Collar(103.0, 22.0, 0.0, longitude=103.0, latitude=22.0),
                total_depth_md=30.0,
                interval_logs={"lithology": [LogInterval(0.0, 30.0, code="C")]},
            )
        ],
    )
    frame = survey_frame(
        [22.0, 22.001, 22.002], [103.0, 103.001, 103.002], ["L", "L", "L"]
    )
    alignment = align_boreholes_to_scene(document, frame, datum="zero")
    traces = scene_borehole_traces(alignment, as_tubes=True)
    assert any(trace.type == "mesh3d" for trace in traces)


def test_borehole_patch_traces_stamps_each_interval_at_true_depth():
    from pycsamt.app._borehole import borehole_patch_traces
    from pycsamt.map.borehole_align import align_boreholes_to_scene
    from pycsamt.map.geometry import survey_frame

    document = PCBHDocument(
        document_id="test:patch",
        created_at="2026-09-04T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:4326"),
        lithologies=[
            VocabularyEntry("A", "Host", color="#f2e4c4"),
            VocabularyEntry("B", "Ore", color="#d59a3f"),
        ],
        boreholes=[
            PCBHBorehole(
                id="BH-1",
                name="Hole",
                kind="mining_exploration",
                status="completed",
                collar=Collar(
                    103.0, 22.0, 0.0, longitude=103.0, latitude=22.0
                ),
                total_depth_md=100.0,
                interval_logs={
                    "lithology": [
                        LogInterval(0.0, 60.0, code="A"),
                        LogInterval(60.0, 100.0, code="B"),
                    ]
                },
            )
        ],
    )
    frame = survey_frame(
        [22.0, 22.001, 22.002], [103.0, 103.001, 103.002], ["L", "L", "L"]
    )
    alignment = align_boreholes_to_scene(document, frame, datum="zero")
    patches = borehole_patch_traces(alignment, half_width=15.0)

    assert len(patches) == 2
    assert all(p.type == "mesh3d" for p in patches)
    colors = {p.color for p in patches}
    assert colors == {"#f2e4c4", "#d59a3f"}
    ore = next(p for p in patches if p.color == "#d59a3f")
    # the ore decal spans its own [from_md, to_md] = [60, 100] m, not the
    # whole hole, and is centred on the collar's along-strike position
    assert min(ore.z) == pytest.approx(-100.0, abs=1e-6)
    assert max(ore.z) == pytest.approx(-60.0, abs=1e-6)
    cx = alignment.placed[0].collar_scene[0]
    assert min(ore.x) == pytest.approx(cx - 15.0, abs=1e-6)
    assert max(ore.x) == pytest.approx(cx + 15.0, abs=1e-6)
