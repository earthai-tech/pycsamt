"""Application-neutral PCBH builder tests."""

from __future__ import annotations

import json

from pycsamt.format.borehole import (
    CSVMappingProfile,
    document_from_builder,
    document_to_builder,
    new_builder_draft,
    pcbh_to_dict,
    preview_csv_mapping,
    validate_builder,
)


def _draft():
    draft = new_builder_draft()
    draft["project"].update(
        {
            "document_id": "builder:test",
            "crs_horizontal": "EPSG:32629",
        }
    )
    draft["boreholes"] = [
        {
            "id": "BH-1",
            "name": "Builder hole",
            "kind": "water",
            "status": "completed",
            "x": 500000,
            "y": 600000,
            "z": 100,
            "total_depth_md": 50,
            "diameter": 0.2,
            "north_reference": "grid",
        }
    ]
    draft["surveys"] = [
        {"borehole_id": "BH-1", "md": 0, "azimuth_deg": 0, "inclination_deg": 0},
        {"borehole_id": "BH-1", "md": 50, "azimuth_deg": 10, "inclination_deg": 15},
    ]
    draft["intervals"] = [
        {
            "borehole_id": "BH-1",
            "family": "lithology",
            "from_md": 0,
            "to_md": 50,
            "code": "SAND",
            "label": "Sand",
            "color": "#D4B483",
        }
    ]
    draft["structures"] = [
        {"borehole_id": "BH-1", "kind": "fracture", "at_md": 25}
    ]
    draft["water"] = [
        {"borehole_id": "BH-1", "at_md": 12, "kind": "water_strike"}
    ]
    draft["construction"] = [
        {"borehole_id": "BH-1", "from_md": 0, "to_md": 20, "kind": "casing"}
    ]
    draft["samples"] = [
        {"borehole_id": "BH-1", "sample_id": "S-1", "from_md": 10, "to_md": 11}
    ]
    draft["assays"] = [
        {"borehole_id": "BH-1", "sample_id": "S-1", "analyte": "Au", "value": 1.2}
    ]
    return draft


def test_builder_constructs_same_canonical_document_and_extensions():
    document = document_from_builder(_draft())
    hole = document.boreholes[0]

    assert document.document_id == "builder:test"
    assert hole.trajectory.method == "survey"
    assert hole.interval_logs["lithology"][0].code == "SAND"
    assert hole.extensions["pcbh:water"][0]["at_md"] == 12
    assert hole.extensions["pcbh:assays"][0]["analyte"] == "Au"


def test_builder_round_trip_is_browser_json_safe():
    document = document_from_builder(_draft())
    restored = document_from_builder(document_to_builder(document))

    assert pcbh_to_dict(restored) == pcbh_to_dict(document)
    json.dumps(document_to_builder(document))


def test_live_validation_routes_issue_to_editor():
    draft = _draft()
    draft["intervals"][0]["to_md"] = 70
    result = validate_builder(draft)

    assert not result.ok
    assert any(item.editor == "intervals" for item in result.diagnostics)
    assert any("total_depth" in item.code or "beyond" in item.code for item in result.diagnostics)


def test_csv_mapping_preview_and_profile_round_trip():
    raw = (
        b"hole,easting,northing,elevation,from,to,lithology,crs\n"
        b"BH-1,1,2,3,0,10,Sand,EPSG:32629\n"
    )
    profile = CSVMappingProfile(
        "mine-columns",
        columns={"borehole.id": "hole"},
    )
    restored = CSVMappingProfile.from_dict(profile.to_dict())
    preview = preview_csv_mapping(raw, profile=restored)

    assert preview["rows"][0]["hole"] == "BH-1"
    assert preview["mapping"]["borehole.id"] == "hole"
    assert preview["mapping"]["collar.x"] == "easting"


def test_empty_builder_reports_actionable_error():
    result = validate_builder(new_builder_draft())
    assert not result.ok
    assert result.diagnostics
