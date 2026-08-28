"""Tests for manifest-driven PCBH relational CSV projects."""

from __future__ import annotations

import yaml
import pytest

from pycsamt.format.borehole import (
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHCSVImportError,
    PCBHDocument,
    RelationalManifest,
    StructureObservation,
    SurveyStation,
    Trajectory,
    VocabularyEntry,
    boreholes_from_csv_directory,
    write_csv_directory,
)


def _document() -> PCBHDocument:
    hole = PCBHBorehole(
        id="BH-1",
        name="BH-1",
        kind="mining_exploration",
        status="completed",
        collar=Collar(100.0, 200.0, 50.0),
        total_depth_md=100.0,
        trajectory=Trajectory(
            method="survey",
            north_reference="grid",
            stations=[
                SurveyStation(0.0, 0.0, 0.0),
                SurveyStation(100.0, 20.0, 10.0),
            ],
        ),
        interval_logs={
            "lithology": [
                LogInterval(0.0, 100.0, code="GRAN", label="Granite")
            ]
        },
        structures=[StructureObservation(kind="fault", at_md=25.0)],
    )
    return PCBHDocument(
        document_id="test:relational",
        created_at="2026-08-28T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:32629"),
        boreholes=[hole],
        lithologies=[VocabularyEntry("GRAN", "Granite")],
        extensions={
            "pcbh:samples": [
                {
                    "sample_id": "S-1",
                    "borehole_id": "BH-1",
                    "from_md": 10.0,
                    "to_md": 12.0,
                    "type": "core",
                }
            ],
            "pcbh:assays": [
                {
                    "sample_id": "S-1",
                    "analyte": "Au",
                    "value": "1.2",
                    "unit": "g/t",
                }
            ],
        },
    )


def test_relational_export_import_preserves_supported_tables(tmp_path):
    root = write_csv_directory(_document(), tmp_path / "project")
    restored, report = boreholes_from_csv_directory(root)

    restored.validate()
    hole = restored.boreholes[0]
    assert len(hole.trajectory.stations) == 2
    assert hole.interval_logs["lithology"][0].label == "Granite"
    assert hole.structures[0].kind == "fault"
    assert restored.extensions["pcbh:samples"][0]["sample_id"] == "S-1"
    assert restored.extensions["pcbh:assays"][0]["analyte"] == "Au"
    assert "collars.csv" in report.source_files
    assert report.ok


def test_missing_optional_tables_are_allowed(tmp_path):
    root = tmp_path / "minimal"
    root.mkdir()
    (root / "collars.csv").write_text(
        "borehole_id,x,y,z,total_depth_md,kind\nB1,1,2,3,10,water\n",
        encoding="utf-8",
    )
    manifest = {
        "version": "0.1.0",
        "crs": {"horizontal": "LOCAL:test"},
        "units": {"depth": "m", "resistivity": "ohm.m"},
        "tables": {"collars": "collars.csv"},
    }
    (root / "import.yaml").write_text(
        yaml.safe_dump(manifest), encoding="utf-8"
    )

    document, _ = boreholes_from_csv_directory(root)
    assert document.boreholes[0].trajectory.method == "vertical"
    assert document.boreholes[0].kind == "water"


def test_broken_join_and_duplicate_sample_fail(tmp_path):
    root = write_csv_directory(_document(), tmp_path / "broken")
    (root / "samples.csv").write_text(
        "sample_id,borehole_id,from_md,to_md\n"
        "S1,MISSING,0,1\nS1,BH-1,1,2\nS1,BH-1,2,3\n",
        encoding="utf-8",
    )
    with pytest.raises(PCBHCSVImportError) as caught:
        boreholes_from_csv_directory(root)
    assert caught.value.report.rows_rejected >= 2


def test_manifest_rejects_traversal_and_unknown_tables(tmp_path):
    bad = {
        "version": "0.1.0",
        "crs": {"horizontal": "LOCAL:test"},
        "tables": {"collars": "../collars.csv", "mystery": "x.csv"},
    }
    path = tmp_path / "import.yaml"
    path.write_text(yaml.safe_dump(bad), encoding="utf-8")
    with pytest.raises(ValueError, match="safe filename|unknown"):
        RelationalManifest.read(path)
