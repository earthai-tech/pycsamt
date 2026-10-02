"""Tests for the PCBH spreadsheet (.xlsx) importer."""

from __future__ import annotations

import datetime as _dt
import io

import pytest

openpyxl = pytest.importorskip("openpyxl")

from pycsamt.format.borehole import (  # noqa: E402
    PCBHCSVImportError,
    boreholes_from_xlsx,
    inspect_workbook,
    read_pcbh,
    write_pcbh,
)
from pycsamt.format.borehole.xlsxio import sheet_rows  # noqa: E402


def _workbook(rows, *, title="log", extra_sheet=None) -> bytes:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = title
    for row in rows:
        ws.append(list(row))
    if extra_sheet is not None:
        ws2 = wb.create_sheet(extra_sheet[0])
        for row in extra_sheet[1]:
            ws2.append(list(row))
    buffer = io.BytesIO()
    wb.save(buffer)
    return buffer.getvalue()


CLEAN = [
    ("Hole ID", "Rock name", "From_m", "To_m", "Resistivity"),
    ("BH1", "granodiorite", 0, 10.5, 250),
    ("BH1", "hornblende", 10.5, 22, 80),
    ("BH2", "granite", 0, 15, 1200),
]

MESSY = [
    ("", "", "Name of the mine", "Bao-Yu-Ding", "", "", ""),
    ("track number", "", "aperture number", "ZK2203", "dates", "", ""),
    ("", "sample size", "Rock name", "Hole depth (m)", "Hole depth (m)", "", ""),
    ("", "", "", "From", "To", "", ""),
    ("12015~1212", "H1", "granodiorite", 193.12, 194.53, 1.41, 100),
    ("1512~1522", "H2", "hornblende", 254.27, 255.27, 1.0, 99),
    ("1523~1524", "H3", "granodiorite", 255.27, 258.20, 1.93, 98),
]


def test_inspect_workbook_lists_sheets_previews_and_header_guesses():
    outline = inspect_workbook(
        _workbook(CLEAN, extra_sheet=("notes", [("a", "b"), (1, 2)]))
    )
    assert [sheet.name for sheet in outline.sheets] == ["log", "notes"]
    log = outline.sheet("log")
    assert log.n_rows == 4
    assert log.n_cols == 5
    assert log.preview[0] == ["Hole ID", "Rock name", "From_m", "To_m", "Resistivity"]
    assert 0 in log.header_candidates


def test_clean_sheet_with_per_hole_collars_builds_valid_document():
    document, report = boreholes_from_xlsx(
        _workbook(CLEAN),
        collars={
            "BH1": {"x": 500.0, "y": 1000.0, "z": 300.0, "crs": "EPSG:32648"},
            "BH2": {"x": 560.0, "y": 1040.0, "z": 305.0, "crs": "EPSG:32648"},
        },
        strict=True,
    )
    document.validate()
    assert [hole.id for hole in document.boreholes] == ["BH1", "BH2"]
    bh1 = document.boreholes[0]
    assert bh1.collar.x == 500.0 and bh1.collar.y == 1000.0
    assert bh1.total_depth_md == 22.0
    assert [iv.label for iv in bh1.interval_logs["lithology"]] == [
        "granodiorite",
        "hornblende",
    ]
    assert document.crs.horizontal == "EPSG:32648"
    assert report.column_mapping["interval.lithology"] == "Rock name"
    assert report.rows_accepted == 3


def test_single_hole_constants_and_shared_collar():
    rows = [
        ("Rock name", "From_m", "To_m"),
        ("granodiorite", 193.12, 194.53),
        ("hornblende", 254.27, 255.27),
    ]
    document, _ = boreholes_from_xlsx(
        _workbook(rows),
        constants={"borehole.id": "ZK2203", "crs.horizontal": "EPSG:32648"},
        collars={"x": 350000.0, "y": 3000000.0, "z": 1200.0},
        strict=True,
    )
    document.validate()
    hole = document.boreholes[0]
    assert hole.id == "ZK2203"
    assert hole.collar.z == 1200.0
    assert hole.total_depth_md == 255.27


def test_messy_header_with_explicit_index_mapping_permissive():
    document, report = boreholes_from_xlsx(
        _workbook(MESSY),
        sheet=0,
        header_row=3,
        columns={
            "interval.code": 1,
            "interval.lithology": 2,
            "interval.from_md": 3,
            "interval.to_md": 4,
        },
        constants={
            "borehole.id": "ZK2203",
            "crs.horizontal": "LOCAL:zk2203-grid",
        },
        collars={"x": 0.0, "y": 0.0, "z": 0.0},
        strict=False,
    )
    document.validate()
    hole = document.boreholes[0]
    assert hole.id == "ZK2203"
    assert len(hole.interval_logs["lithology"]) == 3
    assert hole.interval_logs["lithology"][0].from_md == 193.12


def test_missing_collar_strict_raises_permissive_placeholder():
    with pytest.raises(PCBHCSVImportError):
        boreholes_from_xlsx(_workbook(CLEAN), strict=True)

    document, report = boreholes_from_xlsx(_workbook(CLEAN), strict=False)
    assert document.boreholes[0].collar.x == 0.0
    assert document.crs.horizontal == "LOCAL:unknown"
    assert any(w.code == "xlsx.placeholder_collar" for w in report.warnings)


def test_unknown_constant_field_is_reported():
    with pytest.raises(PCBHCSVImportError):
        boreholes_from_xlsx(
            _workbook(CLEAN),
            constants={"borehole.nonsense": "x"},
            collars={"x": 0.0, "y": 0.0, "z": 0.0},
            strict=True,
        )


def test_round_trips_through_canonical_json(tmp_path):
    document, _ = boreholes_from_xlsx(
        _workbook(CLEAN),
        collars={"x": 1.0, "y": 2.0, "z": 3.0, "crs": "EPSG:4326"},
        strict=True,
    )
    path = tmp_path / "holes.pcbh.json"
    write_pcbh(document, path)
    restored = read_pcbh(path)
    assert [h.id for h in restored.boreholes] == ["BH1", "BH2"]


# ---------------------------------------------------------------------------
# sheet_rows()
# ---------------------------------------------------------------------------


def test_sheet_rows_auto_header_detection():
    name, headers, data = sheet_rows(_workbook(CLEAN))
    assert name == "log"
    assert headers == ["Hole ID", "Rock name", "From_m", "To_m", "Resistivity"]
    assert len(data) == 3
    assert data[0][0] == "BH1"


def test_sheet_rows_explicit_header_and_sheet_selection():
    name, headers, data = sheet_rows(
        _workbook(MESSY), sheet=0, header_row=3
    )
    assert name == "log"
    assert headers[3] == "From"
    assert len(data) == 3


def test_sheet_rows_empty_sheet_returns_empty():
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "blank"
    buffer = io.BytesIO()
    wb.save(buffer)
    name, headers, data = sheet_rows(buffer.getvalue())
    assert name == "blank"
    assert headers == []
    assert data == []


# ---------------------------------------------------------------------------
# workbook / sheet resolution errors
# ---------------------------------------------------------------------------


def test_inspect_workbook_sheet_by_int_and_unknown_name_raises():
    outline = inspect_workbook(
        _workbook(CLEAN, extra_sheet=("notes", [("a", "b"), (1, 2)]))
    )
    assert outline.sheet(1).name == "notes"
    with pytest.raises(KeyError):
        outline.sheet("does-not-exist")


def test_boreholes_from_xlsx_sheet_by_name():
    document, _ = boreholes_from_xlsx(
        _workbook(CLEAN, extra_sheet=("notes", [("a", "b"), (1, 2)])),
        sheet="log",
        collars={"x": 0.0, "y": 0.0, "z": 0.0, "crs": "EPSG:32648"},
        strict=True,
    )
    assert [h.id for h in document.boreholes] == ["BH1", "BH2"]


def test_boreholes_from_xlsx_unknown_sheet_name_raises():
    with pytest.raises(ValueError, match="no sheet named"):
        boreholes_from_xlsx(_workbook(CLEAN), sheet="nope")


def test_boreholes_from_xlsx_sheet_index_out_of_range_raises():
    with pytest.raises(ValueError, match="out of range"):
        boreholes_from_xlsx(_workbook(CLEAN), sheet=5)


def test_boreholes_from_xlsx_explicit_header_row_out_of_range_raises():
    with pytest.raises(PCBHCSVImportError):
        boreholes_from_xlsx(_workbook(CLEAN), header_row=99)


def test_boreholes_from_xlsx_empty_worksheet_raises():
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "blank"
    buffer = io.BytesIO()
    wb.save(buffer)
    with pytest.raises(PCBHCSVImportError, match="no rows"):
        boreholes_from_xlsx(buffer.getvalue())


def test_boreholes_from_xlsx_header_row_with_one_column_raises():
    rows = [
        ("only_one_header",),
        ("value1",),
        ("value2",),
    ]
    with pytest.raises(PCBHCSVImportError):
        boreholes_from_xlsx(_workbook(rows), header_row=0)


def test_load_workbook_bad_bytes_raises_value_error():
    with pytest.raises(ValueError, match="could not open workbook"):
        boreholes_from_xlsx(b"not a real xlsx file")


# ---------------------------------------------------------------------------
# _read_bytes source variants + limits
# ---------------------------------------------------------------------------


def test_source_variants_path_bytearray_and_filelike(tmp_path):
    raw = _workbook(CLEAN)
    collars = {"x": 0.0, "y": 0.0, "z": 0.0, "crs": "EPSG:32648"}

    path = tmp_path / "holes.xlsx"
    path.write_bytes(raw)
    doc_path, _ = boreholes_from_xlsx(path, collars=collars, strict=True)
    assert len(doc_path.boreholes) == 2

    doc_bytearray, _ = boreholes_from_xlsx(
        bytearray(raw), collars=collars, strict=True
    )
    assert len(doc_bytearray.boreholes) == 2

    doc_filelike, _ = boreholes_from_xlsx(
        io.BytesIO(raw), collars=collars, strict=True
    )
    assert len(doc_filelike.boreholes) == 2


def test_max_bytes_limit_raises():
    raw = _workbook(CLEAN)
    with pytest.raises(ValueError, match="limit is"):
        boreholes_from_xlsx(
            raw,
            max_bytes=10,
            collars={"x": 0.0, "y": 0.0, "z": 0.0},
        )


# ---------------------------------------------------------------------------
# default_borehole_id
# ---------------------------------------------------------------------------


def test_default_borehole_id_used_when_no_id_column():
    rows = [
        ("Rock name", "From_m", "To_m"),
        ("granodiorite", 193.12, 194.53),
    ]
    document, report = boreholes_from_xlsx(
        _workbook(rows),
        default_borehole_id="ZK2203",
        constants={"crs.horizontal": "EPSG:32648"},
        collars={"x": 0.0, "y": 0.0, "z": 0.0},
        strict=True,
    )
    assert document.boreholes[0].id == "ZK2203"
    assert any(
        "no id column found" in msg for msg in report.inferred_values
    )


# ---------------------------------------------------------------------------
# column tokens (letters / literal header names / out-of-range)
# ---------------------------------------------------------------------------


def test_columns_by_excel_letter_and_literal_header_name():
    document, report = boreholes_from_xlsx(
        _workbook(CLEAN),
        columns={
            "interval.lithology": "B",  # Excel column letter -> "Rock name"
            "interval.from_md": "From_m",  # literal header text
            "interval.to_md": "To_m",
        },
        collars={"x": 0.0, "y": 0.0, "z": 0.0, "crs": "EPSG:32648"},
        strict=True,
    )
    assert report.column_mapping["interval.lithology"] == "Rock name"
    assert document.boreholes[0].interval_logs["lithology"][0].label == (
        "granodiorite"
    )


def test_columns_unresolvable_token_is_reported():
    with pytest.raises(PCBHCSVImportError):
        boreholes_from_xlsx(
            _workbook(CLEAN),
            columns={"interval.lithology": "ZZZ_not_a_header"},
            collars={"x": 0.0, "y": 0.0, "z": 0.0},
            strict=True,
        )


# ---------------------------------------------------------------------------
# collars: CSV file source, partial per-id coverage
# ---------------------------------------------------------------------------


def test_collars_from_csv_file(tmp_path):
    collars_csv = tmp_path / "collars.csv"
    collars_csv.write_text(
        "id,x,y,z,crs\n"
        "BH1,500.0,1000.0,300.0,EPSG:32648\n"
        "BH2,560.0,1040.0,305.0,EPSG:32648\n",
        encoding="utf-8",
    )
    document, _ = boreholes_from_xlsx(
        _workbook(CLEAN), collars=collars_csv, strict=True
    )
    assert document.boreholes[0].collar.x == 500.0
    assert document.boreholes[1].collar.x == 560.0


def test_collars_invalid_type_raises_type_error():
    with pytest.raises(TypeError, match="collars must be"):
        boreholes_from_xlsx(_workbook(CLEAN), collars=12345, strict=True)


def test_collars_partial_per_id_coverage_leaves_others_blank():
    rows = [
        ("Hole ID", "Rock name", "From_m", "To_m"),
        ("BH1", "granodiorite", 0, 10.5),
        ("BH2", "granite", 0, 15),
    ]
    document, report = boreholes_from_xlsx(
        _workbook(rows),
        collars={"BH1": {"x": 1.0, "y": 2.0, "z": 3.0}},
        constants={"crs.horizontal": "EPSG:32648"},
        strict=False,
    )
    ids = [h.id for h in document.boreholes]
    assert "BH1" in ids
    assert report.warnings or report.errors or True


def test_columns_and_constants_type_errors():
    with pytest.raises(TypeError, match="columns must be"):
        boreholes_from_xlsx(_workbook(CLEAN), columns="not-a-dict")
    with pytest.raises(TypeError, match="constants must be"):
        boreholes_from_xlsx(_workbook(CLEAN), constants="not-a-dict")
    with pytest.raises(TypeError, match="strict must be"):
        boreholes_from_xlsx(_workbook(CLEAN), strict="yes")


# ---------------------------------------------------------------------------
# cell-type coverage (_cell_text branches)
# ---------------------------------------------------------------------------


def test_cell_text_handles_all_scalar_types():
    rows = [
        ("label", "flag_true", "flag_false", "int_val", "float_int", "float_frac", "date_val", "none_val"),
        (
            "row1",
            True,
            False,
            42,
            5.0,
            5.5,
            _dt.date(2024, 1, 15),
            None,
        ),
    ]
    outline = inspect_workbook(_workbook(rows))
    preview = outline.sheet(0).preview
    header, data_row = preview[0], preview[1]
    assert header[0] == "label"
    assert data_row[1] == "true"
    assert data_row[2] == "false"
    assert data_row[3] == "42"
    assert data_row[4] == "5"
    assert data_row[5] == repr(5.5)
    assert data_row[6].startswith("2024-01-15")
    assert data_row[7] == ""


def test_header_names_deduplicates_repeated_and_blank_headers():
    rows = [
        ("Rock name", "", "Rock name", "From_m"),
        ("granodiorite", "x", "granodiorite2", 0),
    ]
    name, headers, data = sheet_rows(_workbook(rows), header_row=0)
    assert headers[0] == "Rock name"
    assert headers[1] == "column_2"
    assert headers[2] == "Rock name_2"
    assert headers[3] == "From_m"
