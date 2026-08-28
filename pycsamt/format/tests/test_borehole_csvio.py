"""Tests for the PCBH combined interval-CSV importer."""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.format.borehole import (
    PCBHCSVImportError,
    boreholes_from_csv,
    read_pcbh,
    write_pcbh,
)

DATA_DIR = Path(__file__).parent / "data" / "borehole"


def test_multi_borehole_template_builds_valid_document_and_report(tmp_path):
    document, report = boreholes_from_csv(
        DATA_DIR / "combined_intervals.csv"
    )

    document.validate()
    assert len(document.boreholes) == 2
    assert report.rows_read == 3
    assert report.rows_accepted == 3
    assert report.rows_rejected == 0
    assert report.ok
    assert len(report.source_sha256) == 64
    assert report.column_mapping["borehole.id"] == "borehole_id"
    assert document.boreholes[0].interval_logs["lithology"][1].to_md == 48.0
    assert document.boreholes[1].interval_logs["lithology"][0].resistivity_ohm_m is None
    assert {entry.name for entry in document.lithologies} == {
        "Laterite",
        "Saprolite",
        "Soil",
    }
    output = tmp_path / "imported.pcbh.json"
    write_pcbh(document, output)
    assert len(read_pcbh(output).boreholes) == 2


def test_explicit_mapping_constants_and_semicolon_detection(tmp_path):
    source = tmp_path / "mapped.csv"
    source.write_text(
        "HOLE;E;N;RL;START;END;ROCK\n"
        "W1;10;20;30;0;15;Clay\n",
        encoding="utf-8",
    )
    document, report = boreholes_from_csv(
        source,
        columns={
            "borehole.id": "HOLE",
            "collar.x": "E",
            "collar.y": "N",
            "collar.z": "RL",
            "interval.from_md": "START",
            "interval.to_md": "END",
            "interval.lithology": "ROCK",
        },
        constants={"crs.horizontal": "LOCAL:test"},
    )

    assert report.delimiter == ";"
    assert document.crs.horizontal == "LOCAL:test"
    assert document.boreholes[0].total_depth_md == 15.0
    assert any("total_depth_md" in item for item in report.inferred_values)


def test_alias_inferences_are_visible_in_report(tmp_path):
    source = tmp_path / "aliases.csv"
    source.write_text(
        "holeid,easting,northing,rl,crs,top,bottom,rock\n"
        "B1,1,2,3,EPSG:32629,0,10,Sand\n",
        encoding="utf-8",
    )
    _, report = boreholes_from_csv(source)

    assert report.warnings
    assert {item.code for item in report.warnings} == {"csv.alias_mapping"}
    assert report.inferred_values


def test_ambiguous_aliases_require_explicit_mapping(tmp_path):
    source = tmp_path / "ambiguous.csv"
    source.write_text(
        "holeid,x,easting,y,z,crs,from,to,lithology\n"
        "B1,1,1,2,3,EPSG:32629,0,10,Sand\n",
        encoding="utf-8",
    )
    with pytest.raises(PCBHCSVImportError) as caught:
        boreholes_from_csv(source)

    assert "csv.mapping_ambiguous" in {
        item.code for item in caught.value.report.errors
    }


def test_repeated_collar_conflict_fails_strict_mode(tmp_path):
    source = tmp_path / "conflict.csv"
    source.write_text(
        "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
        "B1,1,2,3,EPSG:32629,0,10,Sand\n"
        "B1,9,2,3,EPSG:32629,10,20,Clay\n",
        encoding="utf-8",
    )
    with pytest.raises(PCBHCSVImportError) as caught:
        boreholes_from_csv(source)

    report = caught.value.report
    assert report.rows_accepted == 1
    assert report.rows_rejected == 1
    assert report.conflict_resolutions


def test_permissive_mode_returns_partial_document_and_rejections(tmp_path):
    source = tmp_path / "partial.csv"
    source.write_text(
        "borehole_id,x,y,z,crs,from_md,to_md,lithology,resistivity\n"
        "B1,1,2,3,EPSG:32629,0,10,Sand,nan\n"
        "B1,9,2,3,EPSG:32629,10,20,Clay,20\n"
        "B2,4,5,6,EPSG:32629,0,8,nan,30\n",
        encoding="utf-8",
    )
    document, report = boreholes_from_csv(source, strict=False)

    assert [hole.id for hole in document.boreholes] == ["B1"]
    interval = document.boreholes[0].interval_logs["lithology"][0]
    assert interval.resistivity_ohm_m is None
    assert interval.label == "Sand"
    assert report.rows_accepted == 1
    assert report.rows_rejected == 2
    assert len(report.errors) == 2


def test_lithology_codes_are_reused_case_insensitively(tmp_path):
    source = tmp_path / "dictionary.csv"
    source.write_text(
        "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
        "B1,1,2,3,EPSG:32629,0,10,Weathered granite\n"
        "B1,1,2,3,EPSG:32629,10,20,weathered granite\n",
        encoding="utf-8",
    )
    document, _ = boreholes_from_csv(source)

    intervals = document.boreholes[0].interval_logs["lithology"]
    assert intervals[0].code == intervals[1].code
    assert len(document.lithologies) == 1


def test_permissive_rejects_bad_kind_without_adopting_its_crs(tmp_path):
    source = tmp_path / "controlled.csv"
    source.write_text(
        "borehole_id,x,y,z,crs,kind,from_md,to_md,lithology\n"
        "BAD,1,2,3,LOCAL:bad,not_a_kind,0,10,Sand\n"
        "GOOD,4,5,6,EPSG:32629,water,0,10,Clay\n",
        encoding="utf-8",
    )
    document, report = boreholes_from_csv(source, strict=False)

    assert document.crs.horizontal == "EPSG:32629"
    assert [hole.id for hole in document.boreholes] == ["GOOD"]
    assert report.rows_rejected == 1


def test_row_and_byte_resource_limits(tmp_path):
    source = tmp_path / "limited.csv"
    source.write_text(
        "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
        "B1,1,2,3,EPSG:32629,0,10,Sand\n"
        "B1,1,2,3,EPSG:32629,10,20,Clay\n",
        encoding="utf-8",
    )
    with pytest.raises(PCBHCSVImportError) as caught:
        boreholes_from_csv(source, max_rows=1)
    assert caught.value.report.errors[-1].code == "csv.row_limit"

    with pytest.raises(ValueError, match="limit"):
        boreholes_from_csv(source, max_bytes=5)
