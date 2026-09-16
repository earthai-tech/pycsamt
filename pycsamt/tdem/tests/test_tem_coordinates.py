"""Tests for TEM profile/point coordinate readers."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from pycsamt.tdem import (
    TEMCoordinateTable,
    read_tem_coordinates,
    read_temavg_survey,
)
from pycsamt.tdem import coordinates as coords_mod
from pycsamt.tdem.io import (
    read_tem_coordinates as read_tem_coordinates_io,
)

DATA_DIR = Path(__file__).parents[3] / "data" / "TEMAVG" / "JIANGSU"
AVG_FILE = DATA_DIR / "TEM100.AVG"
XLS_FILE = DATA_DIR / "Coordinate of measuring point.xls"


pytestmark = pytest.mark.skipif(
    not AVG_FILE.exists(),
    reason="TEMAVG sample data not available",
)


def _write_coordinate_csv(path: Path) -> None:
    """Create a small coordinate table in the field layout."""
    path.write_text(
        "\n".join(
            [
                "Profile,point,Gauss Coordinate,,Relative coordinate,,H(m),",
                ",,X(m),Y(m),X(m),Y(m),,",
                "100,100,4291789.7679,19510112.9006," "100.0034,100.0151,1102.9537,",
                "100,120,4291809.1,19510120.2,120.0,100.5,1103.1,road",
            ]
        ),
        encoding="utf-8",
    )


def test_read_tem_coordinates_csv(tmp_path):
    """Coordinate CSV files should be parsed by profile and point."""
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)

    coords = read_tem_coordinates(coord_path)
    first = coords.get(100, 100)
    second = coords.get(100, 120)

    assert isinstance(coords, TEMCoordinateTable)
    assert coords.n_points == 2
    assert coords.profiles == [100.0]
    assert coords.points == [100.0, 120.0]
    assert first is not None
    assert first.gauss_x == pytest.approx(4291789.7679)
    assert first.gauss_y == pytest.approx(19510112.9006)
    assert first.x == pytest.approx(100.0034)
    assert first.y == pytest.approx(100.0151)
    assert first.elevation == pytest.approx(1102.9537)
    assert second is not None
    assert second.remark == "road"


def test_read_tem_coordinates_io_wrapper(tmp_path):
    """The public IO module should expose the coordinate reader."""
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)

    coords = read_tem_coordinates_io(coord_path)

    assert coords.n_points == 2
    assert coords.get(100, 100) is not None


def test_survey_records_are_enriched_from_coordinate_file(tmp_path):
    """Explicit coordinates should be attached to matching AVG rows."""
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)

    survey = read_temavg_survey(DATA_DIR, coordinate_file=coord_path)
    coord = survey.coordinate_for(100, 100)
    first = survey.to_records()[0]

    assert survey.coordinates is not None
    assert coord is not None
    assert first["profile"] == pytest.approx(100.0)
    assert first["station"] == pytest.approx(100.0)
    assert first["coord_profile"] == pytest.approx(100.0)
    assert first["coord_point"] == pytest.approx(100.0)
    assert first["x"] == pytest.approx(100.0034)
    assert first["y"] == pytest.approx(100.0151)
    assert first["elevation"] == pytest.approx(1102.9537)


def test_survey_soundings_are_enriched_from_coordinate_file(tmp_path):
    """Generated soundings should carry matching station coordinates."""
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)
    survey = read_temavg_survey(DATA_DIR, coordinate_file=coord_path)

    soundings = survey.to_soundings(stems=["TEM100"])
    first = soundings[0]

    assert first.station_name == "TEM100_100"
    assert first.x == pytest.approx(100.0034)
    assert first.y == pytest.approx(100.0151)
    assert first.elevation == pytest.approx(1102.9537)


# ── TEMCoordinateTable: classmethod / record conversion ──────────────────────


def test_table_read_classmethod_matches_function(tmp_path):
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)

    via_classmethod = TEMCoordinateTable.read(coord_path)
    via_function = read_tem_coordinates(coord_path)

    assert via_classmethod.n_points == via_function.n_points
    assert via_classmethod.coordinates.keys() == via_function.coordinates.keys()


def test_to_records_returns_row_dicts(tmp_path):
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)
    table = read_tem_coordinates(coord_path)

    records = table.to_records()

    assert len(records) == 2
    first = next(r for r in records if r["point"] == 100.0)
    assert first["profile"] == pytest.approx(100.0)
    assert first["gauss_x"] == pytest.approx(4291789.7679)
    assert first["remark"] == ""
    second = next(r for r in records if r["point"] == 120.0)
    assert second["remark"] == "road"


def test_to_dataframe_raises_import_error_without_pandas(monkeypatch, tmp_path):
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)
    table = read_tem_coordinates(coord_path)

    monkeypatch.setitem(sys.modules, "pandas", None)
    with pytest.raises(ImportError, match="requires pandas"):
        table.to_dataframe()


def test_to_dataframe_returns_dataframe(tmp_path):
    pd = pytest.importorskip("pandas")
    coord_path = tmp_path / "coordinates.csv"
    _write_coordinate_csv(coord_path)
    table = read_tem_coordinates(coord_path)

    df = table.to_dataframe()

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2
    assert {"profile", "point", "gauss_x", "gauss_y", "x", "y", "elevation"} <= set(
        df.columns
    )


# ── read_tem_coordinates: file format / error branches ───────────────────────


def test_read_tem_coordinates_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_tem_coordinates(tmp_path / "does-not-exist.csv")


def test_read_tem_coordinates_unsupported_suffix_raises(tmp_path):
    bad = tmp_path / "coordinates.doc"
    bad.write_text("irrelevant", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported TEM coordinate file suffix"):
        read_tem_coordinates(bad)


def test_read_tem_coordinates_tsv(tmp_path):
    coord_path = tmp_path / "coordinates.tsv"
    coord_path.write_text(
        "100\t100\t4291789.7679\t19510112.9006\t100.0034\t100.0151\t1102.9537\t\n"
        "100\t120\t4291809.1\t19510120.2\t120.0\t100.5\t1103.1\troad\n",
        encoding="utf-8",
    )
    table = read_tem_coordinates(coord_path)
    assert table.n_points == 2
    first = table.get(100, 100)
    assert first.gauss_x == pytest.approx(4291789.7679)


def test_read_tem_coordinates_real_xls_sample():
    pytest.importorskip("pandas")
    try:
        table = read_tem_coordinates(XLS_FILE)
    except ImportError:
        pytest.skip("xlrd not available for legacy .xls and no LibreOffice on PATH")
    assert table.n_points > 0


def test_read_tem_coordinates_xlsx(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    pytest.importorskip("pandas")
    xlsx_path = tmp_path / "coordinates.xlsx"
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append([100, 100, 4291789.7679, 19510112.9006, 100.0034, 100.0151, 1102.9537, ""])
    ws.append([100, 120, 4291809.1, 19510120.2, 120.0, 100.5, 1103.1, "road"])
    wb.save(xlsx_path)

    table = read_tem_coordinates(xlsx_path)

    assert table.n_points == 2
    second = table.get(100, 120)
    assert second.remark == "road"


# ── _read_excel_rows: pandas branch selection ─────────────────────────────────


def test_read_excel_rows_xlsx_propagates_pandas_exception(monkeypatch, tmp_path):
    import pandas as pd

    path = tmp_path / "broken.xlsx"
    path.write_text("not a real workbook", encoding="utf-8")

    def _boom(*a, **k):
        raise ValueError("corrupt workbook")

    monkeypatch.setattr(pd, "read_excel", _boom)
    with pytest.raises(ValueError, match="corrupt workbook"):
        coords_mod._read_excel_rows(path)


def test_read_excel_rows_xls_falls_back_on_pandas_exception(monkeypatch, tmp_path):
    import pandas as pd

    path = tmp_path / "legacy.xls"
    path.write_text("not a real workbook", encoding="utf-8")

    def _boom(*a, **k):
        raise ValueError("cannot parse legacy .xls")

    fallback_rows = [["100", "100", "1", "2", "3", "4", "5", ""]]

    monkeypatch.setattr(pd, "read_excel", _boom)
    monkeypatch.setattr(
        coords_mod, "_convert_xls_to_csv_rows", lambda p: fallback_rows
    )
    rows = coords_mod._read_excel_rows(path)
    assert rows == fallback_rows


def test_read_excel_rows_falls_back_when_pandas_unavailable(monkeypatch, tmp_path):
    path = tmp_path / "legacy.xls"
    fallback_rows = [["100", "100", "1", "2", "3", "4", "5", ""]]

    monkeypatch.setitem(sys.modules, "pandas", None)
    monkeypatch.setattr(
        coords_mod, "_convert_xls_to_csv_rows", lambda p: fallback_rows
    )
    rows = coords_mod._read_excel_rows(path)
    assert rows == fallback_rows


# ── _convert_xls_to_csv_rows: LibreOffice fallback ────────────────────────────


def test_convert_xls_no_libreoffice_raises_import_error(monkeypatch, tmp_path):
    monkeypatch.setattr(coords_mod.shutil, "which", lambda _name: None)
    with pytest.raises(ImportError, match="LibreOffice"):
        coords_mod._convert_xls_to_csv_rows(tmp_path / "legacy.xls")


def test_convert_xls_success_reads_generated_csv(monkeypatch, tmp_path):
    src = tmp_path / "legacy.xls"
    src.write_bytes(b"\x00")

    def _fake_run(cmd, check, capture_output, text, env):
        outdir = Path(env["XDG_RUNTIME_DIR"]).parent
        (outdir / f"{src.stem}.csv").write_text(
            "100,100,1.0,2.0,3.0,4.0,5.0,note\n", encoding="utf-8"
        )
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(coords_mod.shutil, "which", lambda _n: "C:/fake/libreoffice")
    monkeypatch.setattr(coords_mod.subprocess, "run", _fake_run)

    rows = coords_mod._convert_xls_to_csv_rows(src)

    assert rows == [["100", "100", "1.0", "2.0", "3.0", "4.0", "5.0", "note"]]


def test_convert_xls_finds_csv_via_glob_when_stem_differs(monkeypatch, tmp_path):
    src = tmp_path / "legacy.xls"
    src.write_bytes(b"\x00")

    def _fake_run(cmd, check, capture_output, text, env):
        outdir = Path(env["XDG_RUNTIME_DIR"]).parent
        (outdir / "unexpected-name.csv").write_text(
            "100,100,1.0,2.0,3.0,4.0,5.0,\n", encoding="utf-8"
        )
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(coords_mod.shutil, "which", lambda _n: "C:/fake/libreoffice")
    monkeypatch.setattr(coords_mod.subprocess, "run", _fake_run)

    rows = coords_mod._convert_xls_to_csv_rows(src)

    assert rows[0][0] == "100"


def test_convert_xls_no_csv_created_raises_file_not_found(monkeypatch, tmp_path):
    src = tmp_path / "legacy.xls"
    src.write_bytes(b"\x00")

    def _fake_run(cmd, check, capture_output, text, env):
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(coords_mod.shutil, "which", lambda _n: "C:/fake/libreoffice")
    monkeypatch.setattr(coords_mod.subprocess, "run", _fake_run)

    with pytest.raises(FileNotFoundError, match="did not create a CSV"):
        coords_mod._convert_xls_to_csv_rows(src)


def test_convert_xls_subprocess_failure_raises_runtime_error(monkeypatch, tmp_path):
    src = tmp_path / "legacy.xls"
    src.write_bytes(b"\x00")

    def _fake_run(cmd, check, capture_output, text, env):
        raise subprocess.CalledProcessError(
            1, cmd, output="", stderr="soffice: fatal error"
        )

    monkeypatch.setattr(coords_mod.shutil, "which", lambda _n: "C:/fake/libreoffice")
    monkeypatch.setattr(coords_mod.subprocess, "run", _fake_run)

    with pytest.raises(RuntimeError, match="soffice: fatal error"):
        coords_mod._convert_xls_to_csv_rows(src)


# ── _parse_coordinate_rows: row-level edge cases ──────────────────────────────


def test_parse_coordinate_rows_skips_short_rows():
    rows = [
        ["100", "100", "1", "2", "3", "4"],  # too short (6 cols)
        ["100", "120", "1", "2", "3", "4", "5", "ok"],
    ]
    coords = coords_mod._parse_coordinate_rows(rows)
    assert list(coords.keys()) == [(100.0, 120.0)]


def test_parse_coordinate_rows_skips_non_numeric_rows():
    rows = [
        ["not-a-number", "100", "1", "2", "3", "4", "5"],
        ["100", "120", "1", "2", "3", "4", "5", "ok"],
    ]
    coords = coords_mod._parse_coordinate_rows(rows)
    assert list(coords.keys()) == [(100.0, 120.0)]


def test_parse_coordinate_rows_all_invalid_raises_value_error():
    rows = [["a", "b", "c", "d", "e", "f", "g"]]
    with pytest.raises(ValueError, match="No coordinate rows found"):
        coords_mod._parse_coordinate_rows(rows)


def test_parse_coordinate_rows_empty_list_raises_value_error():
    with pytest.raises(ValueError, match="No coordinate rows found"):
        coords_mod._parse_coordinate_rows([])
