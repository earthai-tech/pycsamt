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
