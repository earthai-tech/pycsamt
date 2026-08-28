# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Phase 1 tests for the PCBH in-memory and JSON schemas."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pycsamt.format.borehole import (
    PCBH_VERSION,
    Collar,
    CoordinateReferenceSystem,
    LogInterval,
    PCBHBorehole,
    PCBHDocument,
    PCBHValidationError,
    StructureObservation,
    SurveyStation,
    Trajectory,
    ValidationIssue,
    VocabularyEntry,
)

DATA_DIR = Path(__file__).parent / "data" / "borehole"
SCHEMA_PATH = (
    Path(__file__).parents[1]
    / "borehole"
    / "schemas"
    / "pcbh-0.1.schema.json"
)


def _hole(**overrides) -> PCBHBorehole:
    values = {
        "id": "BH-001",
        "name": "BH-001",
        "kind": "water",
        "status": "completed",
        "collar": Collar(x=421550.2, y=1048230.7, z=347.5),
        "total_depth_md": 80.0,
        "trajectory": Trajectory(method="vertical"),
        "interval_logs": {
            "lithology": [
                LogInterval(
                    from_md=0.0,
                    to_md=80.0,
                    code="GRAN",
                    data_nature="observed",
                )
            ]
        },
    }
    values.update(overrides)
    return PCBHBorehole(**values)


def _document(**overrides) -> PCBHDocument:
    values = {
        "document_id": "example:test",
        "created_at": "2026-08-28T12:00:00Z",
        "created_by": "pytest",
        "crs": CoordinateReferenceSystem(horizontal="EPSG:32629"),
        "boreholes": [_hole()],
        "lithologies": [
            VocabularyEntry(code="GRAN", name="Granite", color="#E8A96B")
        ],
    }
    values.update(overrides)
    return PCBHDocument(**values)


class TestValidationIssue:
    def test_rejects_unknown_severity(self):
        with pytest.raises(ValueError, match="severity"):
            ValidationIssue("fatal", "x", "bad")


class TestCoordinateAndCollar:
    def test_valid_coordinates(self):
        crs = CoordinateReferenceSystem(horizontal="LOCAL:mine-grid")
        crs.validate()
        Collar(x=1.0, y=2.0, z=3.0, longitude=-5.0, latitude=9.0).validate()

    def test_lon_lat_must_be_paired(self):
        issues = Collar(x=1.0, y=2.0, z=3.0, longitude=-5.0).collect_issues()
        assert "collar.lonlat_pair" in {issue.code for issue in issues}

    def test_crs_uses_xy_axis_order(self):
        crs = CoordinateReferenceSystem(horizontal="EPSG:4326", axis_order="yx")
        with pytest.raises(PCBHValidationError):
            crs.validate()


class TestTrajectory:
    def test_vertical_has_no_stations(self):
        Trajectory(method="vertical").validate()
        bad = Trajectory(
            method="vertical",
            stations=[SurveyStation(0.0, 0.0, 0.0)],
        )
        assert "trajectory.vertical_stations" in {
            issue.code for issue in bad.collect_issues()
        }

    def test_survey_requires_increasing_stations(self):
        trajectory = Trajectory(
            method="survey",
            north_reference="grid",
            stations=[
                SurveyStation(0.0, 10.0, 5.0),
                SurveyStation(0.0, 12.0, 6.0),
            ],
        )
        assert "trajectory.md_order" in {
            issue.code for issue in trajectory.collect_issues()
        }

    def test_reentry_is_a_warning_not_error(self):
        station = SurveyStation(10.0, 20.0, 100.0)
        issues = station.collect_issues()
        assert [(issue.code, issue.severity) for issue in issues] == [
            ("trajectory.reentry", "warning")
        ]
        station.validate()


class TestIntervalsAndStructures:
    def test_interval_uses_positive_half_open_bounds(self):
        LogInterval(0.0, 10.0, label="Soil").validate()
        with pytest.raises(PCBHValidationError):
            LogInterval(10.0, 10.0, label="Soil").validate()

    def test_interval_needs_code_or_label(self):
        interval = LogInterval(0.0, 10.0)
        assert "interval.identity" in {
            issue.code for issue in interval.collect_issues()
        }

    def test_global_plane_requires_dip_and_direction(self):
        structure = StructureObservation(
            kind="fracture",
            at_md=5.0,
            orientation_representation="global_plane",
            dip_deg=45.0,
        )
        assert "structure.orientation_required" in {
            issue.code for issue in structure.collect_issues()
        }

    def test_structure_location_is_point_xor_interval(self):
        structure = StructureObservation(
            kind="fault",
            at_md=5.0,
            from_md=4.0,
            to_md=6.0,
        )
        assert "structure.location" in {
            issue.code for issue in structure.collect_issues()
        }


class TestBoreholeAndDocument:
    def test_valid_document(self):
        document = _document()
        document.validate()
        assert document.pcbh_version == PCBH_VERSION

    def test_custom_kind_must_be_namespaced(self):
        bad = _hole(kind="blast_hole")
        assert "borehole.kind" in {issue.code for issue in bad.collect_issues()}
        _hole(kind="example:blast_hole").validate()

    def test_depth_members_cannot_exceed_total_depth(self):
        hole = _hole(
            total_depth_md=20.0,
            interval_logs={
                "lithology": [
                    LogInterval(0.0, 30.0, code="GRAN")
                ]
            },
        )
        assert "borehole.interval_beyond_td" in {
            issue.code for issue in hole.collect_issues()
        }

    def test_exclusive_log_family_rejects_overlap(self):
        hole = _hole(
            interval_logs={
                "lithology": [
                    LogInterval(0.0, 12.0, code="GRAN"),
                    LogInterval(10.0, 20.0, code="GRAN"),
                ]
            }
        )
        assert "interval.overlap" in {
            issue.code for issue in hole.collect_issues()
        }

    def test_duplicate_ids_and_unresolved_codes(self):
        document = _document(boreholes=[_hole(), _hole()])
        codes = {issue.code for issue in document.collect_issues()}
        assert "borehole.duplicate_id" in codes

        document = _document(lithologies=[])
        assert "interval.unresolved_code" in {
            issue.code for issue in document.collect_issues()
        }

    def test_case_only_id_collision_is_warning(self):
        second = _hole(id="bh-001", name="second")
        issues = _document(boreholes=[_hole(), second]).collect_issues()
        matching = [issue for issue in issues if issue.code == "borehole.case_collision"]
        assert len(matching) == 1
        assert matching[0].severity == "warning"


class TestJsonSchemaAndFixtures:
    @pytest.fixture(scope="class")
    def validator(self):
        jsonschema = pytest.importorskip("jsonschema")
        schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
        jsonschema.Draft202012Validator.check_schema(schema)
        return jsonschema.Draft202012Validator(
            schema,
            format_checker=jsonschema.FormatChecker(),
        )

    @pytest.mark.parametrize(
        "name",
        ["minimal_vertical.pcbh.json", "multiple_deviated.pcbh.json"],
    )
    def test_valid_fixtures_match_json_schema(self, validator, name):
        document = json.loads((DATA_DIR / name).read_text(encoding="utf-8"))
        assert list(validator.iter_errors(document)) == []

    def test_invalid_crs_fixture_fails_json_schema(self, validator):
        document = json.loads(
            (DATA_DIR / "invalid_crs.pcbh.json").read_text(encoding="utf-8")
        )
        paths = [list(error.absolute_path) for error in validator.iter_errors(document)]
        assert ["crs", "horizontal"] in paths
        assert ["crs", "axis_order"] in paths

    def test_overlap_fixture_is_intentionally_schema_valid(self, validator):
        """Interval overlap is a semantic cross-record rule, not JSON shape."""
        document = json.loads(
            (DATA_DIR / "invalid_overlap.pcbh.json").read_text(encoding="utf-8")
        )
        assert list(validator.iter_errors(document)) == []

