"""Tests for the PCBH spreadsheet (.xlsx) importer."""

from __future__ import annotations

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
