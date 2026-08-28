"""Interoperability tests for PCBH visualization exports."""

from __future__ import annotations

import base64
import json
import struct
import xml.etree.ElementTree as ET

import pytest

from pycsamt.format.borehole import (
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    StructureObservation,
    SurveyStation,
    Trajectory,
    VocabularyEntry,
    write_geojson,
    write_gltf,
    write_vtp,
)


def _document() -> PCBHDocument:
    return PCBHDocument(
        document_id="test:external-exports",
        created_at="2026-08-28T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem(
            "EPSG:4326",
            vertical="EPSG:4979",
        ),
        lithologies=[
            VocabularyEntry("CLAY", "Clay", color="#765432"),
            VocabularyEntry("GRAN", "Granite", color="#A0A0A0"),
        ],
        boreholes=[
            PCBHBorehole(
                id="BH-1",
                name="Export hole",
                kind="mining_exploration",
                status="completed",
                collar=Collar(-4.01, 5.02, 120.0),
                total_depth_md=30.0,
                diameter=0.2,
                trajectory=Trajectory(
                    method="survey",
                    north_reference="grid",
                    stations=[
                        SurveyStation(0.0, 0.0, 0.0),
                        SurveyStation(30.0, 45.0, 20.0),
                    ]
                ),
                interval_logs={
                    "lithology": [
                        LogInterval(0.0, 12.0, code="CLAY"),
                        LogInterval(12.0, 30.0, code="GRAN"),
                    ],
                    "weathering": [
                        LogInterval(0.0, 5.0, label="Weathered")
                    ],
                },
                structures=[
                    StructureObservation(kind="fracture", at_md=15.0)
                ],
            )
        ],
    )


def test_geojson_is_rfc7946_wgs84_with_documented_z(tmp_path):
    path, report = write_geojson(_document(), tmp_path / "holes")
    payload = json.loads(path.read_text(encoding="utf-8"))

    assert payload["type"] == "FeatureCollection"
    assert [item["geometry"]["type"] for item in payload["features"]] == [
        "Point",
        "LineString",
    ]
    assert payload["features"][0]["geometry"]["coordinates"] == [
        -4.01,
        5.02,
        120.0,
    ]
    assert payload["features"][1]["properties"]["z_reference"] == "EPSG:4979"
    assert report.target_crs == "EPSG:4326"
    assert any(loss.code == "geojson_z_semantics" for loss in report.losses)


def test_geojson_rejects_non_rfc7946_target(tmp_path):
    with pytest.raises(ValueError, match="EPSG:4326"):
        write_geojson(
            _document(),
            tmp_path / "bad.geojson",
            target_crs="EPSG:3857",
        )


def test_vtp_parses_as_polydata_with_tubes_and_scalars(tmp_path):
    path, report = write_vtp(_document(), tmp_path / "holes.vtp", sides=6)
    root = ET.parse(path).getroot()
    piece = root.find("./PolyData/Piece")

    assert root.attrib["type"] == "PolyData"
    assert int(piece.attrib["NumberOfPoints"]) > 0
    assert int(piece.attrib["NumberOfPolys"]) > 0
    arrays = {
        item.attrib.get("Name")
        for item in piece.findall("./PointData/DataArray")
    }
    assert {"measured_depth", "material_id", "RGB"} <= arrays
    assert report.features == int(piece.attrib["NumberOfPolys"])
    assert any(
        loss.code == "interval_families_omitted" for loss in report.losses
    )


def test_gltf_embeds_valid_buffer_and_mesh_accessors(tmp_path):
    path, report = write_gltf(_document(), tmp_path / "holes.gltf")
    payload = json.loads(path.read_text(encoding="utf-8"))
    uri = payload["buffers"][0]["uri"]
    binary = base64.b64decode(uri.split(",", 1)[1], validate=True)

    assert payload["asset"]["version"] == "2.0"
    assert len(binary) == payload["buffers"][0]["byteLength"]
    primitive = payload["meshes"][0]["primitives"][0]
    assert primitive["attributes"] == {
        "POSITION": 0,
        "NORMAL": 1,
        "COLOR_0": 2,
    }
    assert payload["accessors"][primitive["indices"]]["count"] > 0
    assert report.format == "gltf"


def test_glb_has_standard_header_and_readable_json_chunk(tmp_path):
    path, report = write_gltf(_document(), tmp_path / "holes.glb")
    raw = path.read_bytes()
    magic, version, total = struct.unpack("<4sII", raw[:12])
    json_length, chunk_type = struct.unpack("<I4s", raw[12:20])
    payload = json.loads(raw[20:20 + json_length].decode("utf-8"))

    assert (magic, version, total) == (b"glTF", 2, len(raw))
    assert chunk_type == b"JSON"
    assert payload["asset"]["version"] == "2.0"
    assert report.format == "glb"
