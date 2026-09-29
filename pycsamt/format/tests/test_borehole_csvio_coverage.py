"""Deep coverage tests for pycsamt.format.borehole.csvio.

Exercises resource-limit / malformed-input errors, the group-consistency
and CRS-conflict paths, lithology-code collision handling, numeric edge
cases, and the document.validate() failure -> PCBHCSVImportError bridge
that the happy-path tests in test_borehole_csvio.py do not reach.
"""

from __future__ import annotations

import math

import pytest

from pycsamt.format.borehole import PCBHCSVImportError, boreholes_from_csv
from pycsamt.format.borehole import csvio, schema


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# --------------------------------------------------------------------- #
# top-level input validation
# --------------------------------------------------------------------- #


class TestInputValidation:
    def test_columns_must_be_dict_or_none(self, tmp_path):
        source = _write(tmp_path, "a.csv", "a,b\n1,2\n")
        with pytest.raises(TypeError, match="columns"):
            boreholes_from_csv(source, columns=["not", "a", "dict"])

    def test_constants_must_be_dict_or_none(self, tmp_path):
        source = _write(tmp_path, "a.csv", "a,b\n1,2\n")
        with pytest.raises(TypeError, match="constants"):
            boreholes_from_csv(source, constants=["nope"])

    def test_strict_must_be_bool(self, tmp_path):
        source = _write(tmp_path, "a.csv", "a,b\n1,2\n")
        with pytest.raises(TypeError, match="strict"):
            boreholes_from_csv(source, strict="yes")

    def test_created_by_must_be_non_empty(self, tmp_path):
        source = _write(tmp_path, "a.csv", "a,b\n1,2\n")
        with pytest.raises(ValueError, match="created_by"):
            boreholes_from_csv(source, created_by="   ")

    def test_max_bytes_type_and_range(self, tmp_path):
        source = _write(tmp_path, "a.csv", "a,b\n1,2\n")
        with pytest.raises(TypeError, match="max_bytes"):
            boreholes_from_csv(source, max_bytes=True)
        with pytest.raises(TypeError, match="max_bytes"):
            boreholes_from_csv(source, max_bytes="10")
        with pytest.raises(ValueError, match="max_bytes"):
            boreholes_from_csv(source, max_bytes=0)

    def test_max_rows_type_and_range(self, tmp_path):
        source = _write(tmp_path, "a.csv", "a,b\n1,2\n")
        with pytest.raises(TypeError, match="max_rows"):
            boreholes_from_csv(source, max_rows=True)
        with pytest.raises(ValueError, match="max_rows"):
            boreholes_from_csv(source, max_rows=0)

    def test_non_utf8_bytes_rejected(self, tmp_path):
        source = tmp_path / "latin1.csv"
        source.write_bytes("borehole_id,x\né,1\n".encode("latin-1"))
        with pytest.raises(ValueError, match="UTF-8"):
            boreholes_from_csv(source)

    def test_unknown_constant_field_is_an_error(self, tmp_path):
        source = _write(
            tmp_path,
            "b.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand\n",
        )
        with pytest.raises(PCBHCSVImportError) as caught:
            boreholes_from_csv(source, constants={"totally.unknown": "x"})
        assert "csv.constant_field" in {
            i.code for i in caught.value.report.errors
        }

    def test_empty_csv_file_raises(self, tmp_path):
        source = _write(tmp_path, "empty.csv", "")
        with pytest.raises(PCBHCSVImportError) as caught:
            boreholes_from_csv(source)
        assert caught.value.report.errors[0].code == "csv.empty"

    def test_header_too_short_or_blank_cell(self, tmp_path):
        source = _write(tmp_path, "short_header.csv", "onlyone\n1\n")
        with pytest.raises(PCBHCSVImportError) as caught:
            boreholes_from_csv(source)
        codes = {i.code for i in caught.value.report.errors}
        assert "csv.header" in codes

    def test_explicit_delimiter_must_be_known(self, tmp_path):
        source = _write(tmp_path, "d.csv", "a,b\n1,2\n")
        with pytest.raises(ValueError, match="delimiter must be one of"):
            boreholes_from_csv(source, delimiter="~")

    def test_delimiter_sniff_failure_raises(self, tmp_path):
        # A single-column, single-character body defeats csv.Sniffer.
        source = _write(tmp_path, "nodelim.csv", "x\n")
        with pytest.raises(ValueError, match="delimiter"):
            boreholes_from_csv(source)


# --------------------------------------------------------------------- #
# row-level skip / reject paths
# --------------------------------------------------------------------- #


class TestRowHandling:
    def test_blank_rows_are_skipped_not_rejected(self, tmp_path):
        source = _write(
            tmp_path,
            "blank.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand\n"
            ",,,,,,,\n",
        )
        document, report = boreholes_from_csv(source)
        assert report.rows_skipped == 1
        assert report.rows_accepted == 1
        assert len(document.boreholes) == 1

    def test_row_width_mismatch_is_rejected(self, tmp_path):
        source = _write(
            tmp_path,
            "width.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,10\n",  # missing lithology cell
        )
        document, report = boreholes_from_csv(source, strict=False)
        assert report.rows_rejected == 1
        assert any(i.code == "csv.row_width" for i in report.errors)

    def test_crs_conflict_between_boreholes_is_rejected(self, tmp_path):
        source = _write(
            tmp_path,
            "crs_conflict.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand\n"
            "B2,4,5,6,EPSG:4326,0,10,Clay\n",
        )
        document, report = boreholes_from_csv(source, strict=False)
        assert [h.id for h in document.boreholes] == ["B1"]
        assert any(i.code == "csv.crs_conflict" for i in report.errors)

    def test_interval_overlap_within_borehole_is_rejected(self, tmp_path):
        source = _write(
            tmp_path,
            "overlap.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand\n"
            "B1,1,2,3,EPSG:32629,5,15,Clay\n",
        )
        document, report = boreholes_from_csv(source, strict=False)
        intervals = document.boreholes[0].interval_logs["lithology"]
        assert len(intervals) == 1
        assert any(i.code == "csv.interval_overlap" for i in report.errors)

    def test_interval_beyond_total_depth_is_rejected(self, tmp_path):
        source = _write(
            tmp_path,
            "beyond_td.csv",
            "borehole_id,x,y,z,crs,total_depth_md,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,10,0,5,Sand\n"
            "B1,1,2,3,EPSG:32629,10,5,15,Clay\n",
        )
        document, report = boreholes_from_csv(source, strict=False)
        intervals = document.boreholes[0].interval_logs["lithology"]
        assert len(intervals) == 1
        assert any(i.code == "csv.interval_beyond_td" for i in report.errors)

    def test_all_rows_rejected_raises_no_valid_boreholes(self, tmp_path):
        source = _write(
            tmp_path,
            "allbad.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,10,5,Sand\n",  # to < from -> rejected
        )
        with pytest.raises(PCBHCSVImportError) as caught:
            boreholes_from_csv(source, strict=False)
        assert caught.value.report.errors[-1].code == "csv.no_valid_rows"

    def test_lithology_code_collision_from_explicit_codes(self, tmp_path):
        source = _write(
            tmp_path,
            "code_conflict.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology,code\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand,ROCK1\n"
            "B1,1,2,3,EPSG:32629,10,20,Clay,ROCK1\n",
        )
        document, report = boreholes_from_csv(source, strict=False)
        intervals = document.boreholes[0].interval_logs["lithology"]
        assert len(intervals) == 1
        assert any(
            i.code == "csv.lithology_code_conflict" for i in report.errors
        )

    def test_lithology_code_auto_generation_collision_gets_suffixed(
        self, tmp_path
    ):
        # "Clay!" and "Clay?" both sanitize to the code "CLAY"; the second,
        # different label must be suffixed rather than silently merged.
        source = _write(
            tmp_path,
            "auto_code_collision.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,10,Clay!\n"
            "B1,1,2,3,EPSG:32629,10,20,Clay?\n",
        )
        document, report = boreholes_from_csv(source)
        intervals = document.boreholes[0].interval_logs["lithology"]
        assert intervals[0].code != intervals[1].code
        assert intervals[1].code.startswith("CLAY")

    def test_units_other_than_metres_ohm_metres_rejected(self, tmp_path):
        source = _write(
            tmp_path,
            "units.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand\n",
        )
        with pytest.raises(PCBHCSVImportError) as caught:
            boreholes_from_csv(
                source, constants={"units.depth": "ft"}
            )
        assert any(
            i.code == "csv.units_unsupported" for i in caught.value.report.errors
        )

    def test_document_validation_failure_is_wrapped(self, tmp_path, monkeypatch):
        source = _write(
            tmp_path,
            "valid.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand\n",
        )

        def _boom(self):
            raise schema.PCBHValidationError(
                [
                    schema.ValidationIssue(
                        code="document.fabricated",
                        message="forced failure",
                        path="$.fabricated",
                        severity="error",
                    )
                ]
            )

        monkeypatch.setattr(schema.PCBHDocument, "validate", _boom)
        with pytest.raises(PCBHCSVImportError) as caught:
            boreholes_from_csv(source)
        assert any(
            i.code == "csv.document.document.fabricated"
            for i in caught.value.report.errors
        )


# --------------------------------------------------------------------- #
# _parse_row field-level errors
# --------------------------------------------------------------------- #


class TestParseRowFieldErrors:
    def test_missing_and_blank_hole_id(self, tmp_path):
        source = _write(
            tmp_path,
            "id.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            ",1,2,3,EPSG:32629,0,10,Sand\n"
            "  ,1,2,3,EPSG:32629,0,10,Sand\n",
        )
        _, report = boreholes_from_csv(source, strict=False)
        codes = [i.code for i in report.errors]
        assert codes.count("csv.missing_id") == 2

    def test_missing_and_blank_lithology(self, tmp_path):
        source = _write(
            tmp_path,
            "lith.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,10,\n"
            "B2,1,2,3,EPSG:32629,0,10,  \n",
        )
        _, report = boreholes_from_csv(source, strict=False)
        codes = [i.code for i in report.errors]
        assert codes.count("csv.missing_lithology") == 2

    def test_missing_and_blank_crs(self, tmp_path):
        source = _write(
            tmp_path,
            "crs.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,,0,10,Sand\n"
            "B2,1,2,3,  ,0,10,Clay\n",
        )
        _, report = boreholes_from_csv(source, strict=False)
        codes = [i.code for i in report.errors]
        assert codes.count("csv.missing_crs") == 2

    def test_interval_bounds_invalid(self, tmp_path):
        source = _write(
            tmp_path,
            "bounds.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,-1,10,Sand\n"
            "B2,1,2,3,EPSG:32629,10,10,Clay\n",
        )
        _, report = boreholes_from_csv(source, strict=False)
        assert sum(
            i.code == "csv.interval_bounds" for i in report.errors
        ) == 2

    def test_resistivity_must_be_positive(self, tmp_path):
        source = _write(
            tmp_path,
            "resist.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology,resistivity\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand,0\n",
        )
        _, report = boreholes_from_csv(source, strict=False)
        assert any(i.code == "csv.resistivity" for i in report.errors)

    def test_total_depth_must_be_positive(self, tmp_path):
        source = _write(
            tmp_path,
            "td.csv",
            "borehole_id,x,y,z,crs,total_depth_md,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,0,10,Sand\n",
        )
        _, report = boreholes_from_csv(source, strict=False)
        assert any(i.code == "csv.total_depth" for i in report.errors)

    def test_unsupported_borehole_kind_and_status(self, tmp_path):
        source = _write(
            tmp_path,
            "kindstatus.csv",
            "borehole_id,x,y,z,crs,kind,status,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,not_a_kind,unknown,0,10,Sand\n"
            "B2,1,2,3,EPSG:32629,water,not_a_status,0,10,Clay\n",
        )
        _, report = boreholes_from_csv(source, strict=False)
        codes = {i.code for i in report.errors}
        assert "csv.borehole_kind" in codes
        assert "csv.borehole_status" in codes

    def test_namespaced_kind_and_status_are_accepted(self, tmp_path):
        source = _write(
            tmp_path,
            "namespaced.csv",
            "borehole_id,x,y,z,crs,kind,status,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,vendor:custom_kind,vendor:custom_status,"
            "0,10,Sand\n",
        )
        document, report = boreholes_from_csv(source)
        assert report.ok
        assert document.boreholes[0].kind == "vendor:custom_kind"
        assert document.boreholes[0].status == "vendor:custom_status"

    def test_invalid_data_nature_triggers_interval_collect_issues(
        self, tmp_path
    ):
        source = _write(
            tmp_path,
            "datanature.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology,data_nature\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand,not_a_nature\n",
        )
        with pytest.raises(PCBHCSVImportError) as caught:
            boreholes_from_csv(source)
        assert any(
            i.code == "csv.interval.data_nature"
            for i in caught.value.report.errors
        )

    def test_boolean_numeric_value_is_invalid_number(self, tmp_path):
        # A stray boolean-looking cell for a numeric column must not be
        # silently coerced to 0/1; it should register as invalid.
        source = _write(
            tmp_path,
            "bool_numeric.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology,resistivity\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand,True\n",
        )
        _, report = boreholes_from_csv(
            source,
            constants={"interval.resistivity_ohm_m": True},
            strict=False,
        )
        assert any(i.code == "csv.invalid_number" for i in report.errors)

    def test_non_finite_number_is_invalid(self, tmp_path):
        source = _write(
            tmp_path,
            "inf.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,inf,2,3,EPSG:32629,0,10,Sand\n",
        )
        _, report = boreholes_from_csv(source, strict=False)
        assert any(i.code == "csv.invalid_number" for i in report.errors)


# --------------------------------------------------------------------- #
# group consistency
# --------------------------------------------------------------------- #


class TestGroupConsistency:
    def test_total_depth_supplied_on_a_later_row_is_adopted(self, tmp_path):
        source = _write(
            tmp_path,
            "later_td.csv",
            "borehole_id,x,y,z,crs,total_depth_md,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,,0,10,Sand\n"
            "B1,1,2,3,EPSG:32629,20,10,20,Clay\n",
        )
        document, report = boreholes_from_csv(source)
        assert report.ok
        assert document.boreholes[0].total_depth_md == 20.0

    def test_total_depth_missing_on_every_row_is_inferred(self, tmp_path):
        source = _write(
            tmp_path,
            "no_td.csv",
            "borehole_id,x,y,z,crs,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,0,10,Sand\n"
            "B1,1,2,3,EPSG:32629,10,25,Clay\n",
        )
        document, report = boreholes_from_csv(source)
        assert document.boreholes[0].total_depth_md == 25.0
        assert any(
            "total_depth_md" in item for item in report.inferred_values
        )

    def test_collar_conflict_reports_which_field(self, tmp_path):
        source = _write(
            tmp_path,
            "namefield.csv",
            "borehole_id,x,y,z,crs,name,from_md,to_md,lithology\n"
            "B1,1,2,3,EPSG:32629,First,0,10,Sand\n"
            "B1,1,2,3,EPSG:32629,Second,10,20,Clay\n",
        )
        with pytest.raises(PCBHCSVImportError) as caught:
            boreholes_from_csv(source)
        assert any(
            "borehole.name" in i.message for i in caught.value.report.errors
        )


# --------------------------------------------------------------------- #
# _missing / _controlled helpers directly
# --------------------------------------------------------------------- #


class TestHelpers:
    @pytest.mark.parametrize(
        "value",
        [None, "", "  ", "na", "N/A", "NaN", "none", "NULL"],
    )
    def test_missing_recognizes_sentinels(self, value):
        assert csvio._missing(value) is True

    def test_missing_false_for_real_values(self):
        assert csvio._missing("Sand") is False
        assert csvio._missing(0) is False

    def test_controlled_accepts_namespaced_values(self):
        assert csvio._controlled("water", schema.BOREHOLE_KINDS) is True
        assert csvio._controlled("vendor:x", schema.BOREHOLE_KINDS) is True
        assert csvio._controlled("vendor:", schema.BOREHOLE_KINDS) is False
        assert csvio._controlled("nonsense", schema.BOREHOLE_KINDS) is False

    def test_number_helper_directly(self):
        report = csvio.ImportReport(
            source="x", source_sha256="0" * 64, delimiter=",", strict=True
        )
        assert csvio._number(None, "x", 1, "col", report, required=False) is None
        assert report.errors == ()

        report2 = csvio.ImportReport(
            source="x", source_sha256="0" * 64, delimiter=",", strict=True
        )
        assert csvio._number(None, "x", 1, "col", report2, required=True) is None
        assert len(report2.errors) == 1

        report3 = csvio.ImportReport(
            source="x", source_sha256="0" * 64, delimiter=",", strict=True
        )
        assert (
            csvio._number(math.inf, "x", 1, "col", report3, required=False)
            is None
        )
        assert report3.errors[0].code == "csv.invalid_number"

        report4 = csvio.ImportReport(
            source="x", source_sha256="0" * 64, delimiter=",", strict=True
        )
        assert csvio._number("3.5", "x", 1, "col", report4, required=False) == 3.5
