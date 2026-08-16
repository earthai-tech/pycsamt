# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.interp.export — Oasis Montaj XYZ, LAS, CSV, VTK,
Surfer grid/XYZ."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.interp import export
from pycsamt.interp._base import ResistivityModel
from pycsamt.geology.lithology import Layer, StratigraphicLog


def _log(name="S1", x=100.0, with_nan=True):
    z = np.array([5.0, 15.0, 25.0, 35.0])
    rho = np.log10(np.array([10.0, 150.0, np.nan, 5000.0]))
    if not with_nan:
        rho = np.log10(np.array([10.0, 150.0, 800.0, 5000.0]))
    layers = [
        Layer(top=0.0, bottom=20.0, rho_log10=np.log10(10.0), lithology="Clay"),
        Layer(
            top=20.0,
            bottom=40.0,
            rho_log10=np.log10(150.0),
            lithology="Sand (wet)",
        ),
    ]
    return StratigraphicLog(
        station_name=name,
        station_x=x,
        z_centers=z,
        rho_log10=rho,
        layers=layers,
    )


def _model():
    rho_2d = np.log10(np.array([[10.0, 150.0], [800.0, np.nan], [5000.0, 200.0]]))
    return ResistivityModel.from_array(
        rho_2d,
        x_centers=np.array([0.0, 100.0]),
        z_centers=np.array([5.0, 15.0, 25.0]),
        method="test",
    )


# ─────────────────────────────────────────────────────────────────────────────
# to_oasis_montaj_xyz
# ─────────────────────────────────────────────────────────────────────────────


def test_to_oasis_montaj_xyz_default(tmp_path):
    out = export.to_oasis_montaj_xyz([_log()], tmp_path / "profile.xyz")
    assert out.exists()
    text = out.read_text()
    assert "/ pycsamt.interp" in text
    assert "/ Line S1" in text
    assert "X" in text and "RESD" in text
    # nan row (z=25.0) must be skipped
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("/")]
    assert len(lines) == 3  # 4 depths - 1 nan


def test_to_oasis_montaj_xyz_z_uses_negative_depth_without_elevation(
    tmp_path,
):
    out = export.to_oasis_montaj_xyz([_log()], tmp_path / "p.xyz")
    lines = [ln for ln in out.read_text().splitlines() if ln and not ln.startswith("/")]
    first_z = float(lines[0].split()[2])
    assert first_z == pytest.approx(-5.0)


def test_to_oasis_montaj_xyz_with_elevation(tmp_path):
    out = export.to_oasis_montaj_xyz(
        [_log()], tmp_path / "p.xyz", elevation=np.array([100.0])
    )
    lines = [ln for ln in out.read_text().splitlines() if ln and not ln.startswith("/")]
    first_z = float(lines[0].split()[2])
    assert first_z == pytest.approx(100.0 - 5.0)


def test_to_oasis_montaj_xyz_linear_rho(tmp_path):
    out = export.to_oasis_montaj_xyz([_log()], tmp_path / "p.xyz", log_rho=False)
    lines = [ln for ln in out.read_text().splitlines() if ln and not ln.startswith("/")]
    rho_val = float(lines[0].split()[3])
    assert rho_val == pytest.approx(10.0)


def test_to_oasis_montaj_xyz_custom_channels_header(tmp_path):
    out = export.to_oasis_montaj_xyz(
        [_log()], tmp_path / "p.xyz", channels=["A", "B", "C", "D", "E"]
    )
    text = out.read_text()
    assert "/ A  B  C  D  E" in text


def test_to_oasis_montaj_xyz_lithology_lookup_and_space_replace(tmp_path):
    out = export.to_oasis_montaj_xyz([_log()], tmp_path / "p.xyz")
    lines = [ln for ln in out.read_text().splitlines() if ln and not ln.startswith("/")]
    # z=5.0 and z=15.0 fall in the Clay layer [0,20); the nan row (z=25)
    # is skipped, so the third written row is z=35.0 -> Sand (wet) [20,40)
    lith_first = lines[0].split()[-1]
    assert lith_first == "Clay"
    lith_third = lines[2].split()[-1]
    assert lith_third == "Sand_(wet)"


def test_to_oasis_montaj_xyz_depth_outside_all_layers_blank_lithology(
    tmp_path,
):
    log = _log()
    log.layers = []  # no layers at all -> lithology lookup always misses
    out = export.to_oasis_montaj_xyz([log], tmp_path / "p.xyz")
    lines = [ln for ln in out.read_text().splitlines() if ln and not ln.startswith("/")]
    # trailing lithology field should be empty (line ends after rho value)
    parts = lines[0].split()
    assert len(parts) == 4  # x, y, z, rho — no lithology token


def test_to_oasis_montaj_xyz_multiple_logs(tmp_path):
    out = export.to_oasis_montaj_xyz(
        [_log("S1", 0.0), _log("S2", 100.0)], tmp_path / "p.xyz"
    )
    text = out.read_text()
    assert "/ Line S1" in text
    assert "/ Line S2" in text


def test_to_oasis_montaj_xyz_creates_parent_dirs(tmp_path):
    out = export.to_oasis_montaj_xyz([_log()], tmp_path / "nested" / "dir" / "p.xyz")
    assert out.exists()


# ─────────────────────────────────────────────────────────────────────────────
# to_las
# ─────────────────────────────────────────────────────────────────────────────


def test_to_las_basic(tmp_path):
    out = export.to_las(_log(), tmp_path / "S1.las")
    text = out.read_text()
    assert "~VERSION INFORMATION" in text
    assert "~WELL INFORMATION" in text
    assert "~CURVE INFORMATION" in text
    assert "~A  DEPT" in text
    assert " WELL.           S1:  WELL" in text
    assert "LOG10.OHMM" in text


def test_to_las_custom_well_name_and_company(tmp_path):
    out = export.to_las(
        _log(), tmp_path / "S1.las", well_name="CustomWell", company="ACME"
    )
    text = out.read_text()
    assert "CustomWell" in text
    assert "ACME" in text


def test_to_las_linear_rho_units(tmp_path):
    out = export.to_las(_log(), tmp_path / "S1.las", log_rho=False)
    text = out.read_text()
    assert " RESD.OHMM  " in text
    assert "LOG10.OHMM" not in text


def test_to_las_null_value_written_for_nan(tmp_path):
    out = export.to_las(_log(), tmp_path / "S1.las", null_value=-999.25)
    text = out.read_text()
    assert "-999.25000" in text or "-999.25" in text


def test_to_las_single_depth_step_default(tmp_path):
    z = np.array([10.0])
    rho = np.array([np.log10(100.0)])
    log = StratigraphicLog(
        station_name="S1", station_x=0.0, z_centers=z, rho_log10=rho, layers=[]
    )
    out = export.to_las(log, tmp_path / "single.las")
    text = out.read_text()
    assert " STEP.M          1.0000:  STEP" in text


def test_to_las_lithology_code_hash(tmp_path):
    out = export.to_las(_log(), tmp_path / "S1.las")
    text = out.read_text()
    data_lines = [
        ln for ln in text.splitlines() if ln.strip() and ln.strip()[0].isdigit()
    ]
    assert len(data_lines) >= 1


def test_to_las_creates_parent_dirs(tmp_path):
    out = export.to_las(_log(), tmp_path / "nested" / "S1.las")
    assert out.exists()


# ─────────────────────────────────────────────────────────────────────────────
# to_csv
# ─────────────────────────────────────────────────────────────────────────────


def test_to_csv_basic(tmp_path):
    out = export.to_csv([_log()], tmp_path / "all.csv")
    text = out.read_text()
    lines = text.strip().splitlines()
    assert lines[0] == "station,x_m,depth_m,rho_log10,rho_ohm_m,lithology"
    # 4 depths - 1 nan row skipped = 3 data rows
    assert len(lines) == 1 + 3


def test_to_csv_lithology_column_populated(tmp_path):
    out = export.to_csv([_log()], tmp_path / "all.csv")
    text = out.read_text()
    assert "Clay" in text
    assert "Sand (wet)" in text  # CSV keeps original spacing, unlike XYZ


def test_to_csv_lithology_blank_when_no_layer_match(tmp_path):
    log = _log(with_nan=False)
    log.layers = []  # no layers at all -> every depth falls through unmatched
    out = export.to_csv([log], tmp_path / "all.csv")
    text = out.read_text()
    lines = text.strip().splitlines()[1:]
    assert all(ln.endswith(",") for ln in lines)  # blank trailing lithology


def test_to_csv_multiple_logs(tmp_path):
    out = export.to_csv(
        [_log("S1", 0.0, with_nan=False), _log("S2", 100.0, with_nan=False)],
        tmp_path / "all.csv",
    )
    text = out.read_text()
    lines = text.strip().splitlines()
    assert len(lines) == 1 + 4 + 4  # header + 4 rows each


def test_to_csv_creates_parent_dirs(tmp_path):
    out = export.to_csv([_log()], tmp_path / "nested" / "dir" / "all.csv")
    assert out.exists()


# ─────────────────────────────────────────────────────────────────────────────
# to_vtk
# ─────────────────────────────────────────────────────────────────────────────


def test_to_vtk_basic(tmp_path):
    out = export.to_vtk(_model(), tmp_path / "m.vtk")
    text = out.read_text()
    assert "# vtk DataFile Version 3.0" in text
    assert "DATASET RECTILINEAR_GRID" in text
    assert "DIMENSIONS 2 3 1" in text
    assert "SCALARS log10_rho float 1" in text


def test_to_vtk_linear_rho(tmp_path):
    out = export.to_vtk(_model(), tmp_path / "m.vtk", log_rho=False)
    text = out.read_text()
    assert "10.00000" in text  # 10**log10(10) == 10


def test_to_vtk_custom_field_name(tmp_path):
    out = export.to_vtk(_model(), tmp_path / "m.vtk", field_name="rho_ohm")
    text = out.read_text()
    assert "SCALARS rho_ohm float 1" in text


def test_to_vtk_nan_replaced_with_sentinel(tmp_path):
    out = export.to_vtk(_model(), tmp_path / "m.vtk")
    text = out.read_text()
    assert "-9999.00000" in text


def test_to_vtk_creates_parent_dirs(tmp_path):
    out = export.to_vtk(_model(), tmp_path / "nested" / "m.vtk")
    assert out.exists()


def test_to_vtk_point_data_count(tmp_path):
    out = export.to_vtk(_model(), tmp_path / "m.vtk")
    text = out.read_text()
    assert "POINT_DATA 6" in text  # nx=2, nz=3


# ─────────────────────────────────────────────────────────────────────────────
# to_surfer_grid
# ─────────────────────────────────────────────────────────────────────────────


def _uniform_model():
    """A model already on a perfectly uniform grid, so resampling is a
    lossless identity transform -- lets tests assert exact values."""
    rho_2d = np.array([[2.0, 2.1, 2.2], [2.5, 2.6, 2.7]])
    return ResistivityModel.from_array(
        rho_2d,
        x_centers=np.array([0.0, 100.0, 200.0]),
        z_centers=np.array([10.0, 50.0]),
        method="test",
    )


class _FakeOccamMesh:
    """Minimal native-Occam2D-shaped mesh: 4 uniform 50 m cells, no
    padding, so x_centers land exactly on 25/75/125/175."""

    def __init__(self):
        self.x_nodes = np.array([0.0, 50.0, 100.0, 150.0, 200.0])
        self.z_nodes = np.array([0.0, 20.0, 60.0])


class _FakeOccamData:
    def __init__(self):
        self.offsets = np.array([25.0, 75.0, 125.0, 175.0])
        self.sites = ["S0", "S1", "S2", "S3"]


class _FakeOccamResult:
    """Duck-typed native pycsamt.models.occam2d.results.InversionResult,
    the same pattern pycsamt/interp/tests/test_base.py uses -- proves
    to_surfer_grid/to_surfer_xyz accept a *native* Occam2D result
    directly, not only an already-built ResistivityModel."""

    def __init__(self):
        self.mesh = _FakeOccamMesh()
        self.data = _FakeOccamData()
        self.rho_2d = np.array([[2.0, 2.1, 2.2, 2.3], [2.5, 2.6, 2.7, 2.8]])
        self.final_rms = 1.5


def test_to_surfer_grid_accepts_native_occam2d_result(tmp_path):
    out = export.to_surfer_grid(_FakeOccamResult(), tmp_path / "m.grd")
    assert out.exists()
    assert out.read_text().splitlines()[0] == "DSAA"


def test_to_surfer_xyz_accepts_native_occam2d_result(tmp_path):
    out = export.to_surfer_xyz(_FakeOccamResult(), tmp_path / "m.dat")
    lines = out.read_text().splitlines()
    assert lines[0] == "X\tY\tZ"
    assert len(lines) == 1 + 4 * 2  # header + 4 stations x 2 depths


def test_to_surfer_grid_dsaa_header(tmp_path):
    out = export.to_surfer_grid(_uniform_model(), tmp_path / "m.grd", nx=3, ny=2)
    lines = out.read_text().splitlines()
    assert lines[0] == "DSAA"
    assert lines[1] == "3 2"
    assert lines[2] == "0 200"
    assert lines[3] == "-50 -10"


def test_to_surfer_grid_values_match_source_on_uniform_grid(tmp_path):
    out = export.to_surfer_grid(_uniform_model(), tmp_path / "m.grd", nx=3, ny=2)
    rows = out.read_text().splitlines()[5:]
    deep_row = [float(v) for v in rows[0].split()]   # y_min = -50 (z=50)
    shallow_row = [float(v) for v in rows[1].split()]  # y_max = -10 (z=10)
    np.testing.assert_allclose(deep_row, [2.5, 2.6, 2.7])
    np.testing.assert_allclose(shallow_row, [2.0, 2.1, 2.2])


def test_to_surfer_grid_linear_rho(tmp_path):
    out = export.to_surfer_grid(
        _uniform_model(), tmp_path / "m.grd", nx=3, ny=2, log_rho=False,
    )
    rows = out.read_text().splitlines()[5:]
    shallow_row = [float(v) for v in rows[1].split()]
    # rtol matches the file's own %.6g write precision, not full float64
    np.testing.assert_allclose(
        shallow_row, 10.0 ** np.array([2.0, 2.1, 2.2]), rtol=1e-5
    )


def test_to_surfer_grid_nan_source_cell_blanked_not_literal_nan(tmp_path):
    out = export.to_surfer_grid(_model(), tmp_path / "m.grd")
    text = out.read_text()
    assert "nan" not in text.lower()


def test_to_surfer_grid_accepts_raw_array_triple(tmp_path):
    x = np.array([0.0, 100.0, 200.0])
    z = np.array([10.0, 50.0])
    rho = np.array([[2.0, 2.1, 2.2], [2.5, 2.6, 2.7]])
    out = export.to_surfer_grid((x, z, rho), tmp_path / "m.grd", nx=3, ny=2)
    assert out.exists()
    assert out.read_text().splitlines()[0] == "DSAA"


def test_to_surfer_grid_with_elevation_shifts_y_range(tmp_path):
    m = _uniform_model()
    flat = export.to_surfer_grid(m, tmp_path / "flat.grd", nx=3, ny=2)
    draped = export.to_surfer_grid(
        m, tmp_path / "draped.grd", nx=3, ny=2,
        elevation=np.array([100.0, 100.0, 100.0]),
        chainage=m.x_centers,
    )
    flat_ylim = flat.read_text().splitlines()[3]
    draped_ylim = draped.read_text().splitlines()[3]
    assert flat_ylim != draped_ylim
    # elevation=100 uniformly shifts Y by +100 relative to the flat case
    y0_flat, y1_flat = (float(v) for v in flat_ylim.split())
    y0_draped, y1_draped = (float(v) for v in draped_ylim.split())
    assert y0_draped == pytest.approx(y0_flat + 100.0)
    assert y1_draped == pytest.approx(y1_flat + 100.0)


def test_to_surfer_grid_creates_parent_dirs(tmp_path):
    out = export.to_surfer_grid(_uniform_model(), tmp_path / "nested" / "m.grd")
    assert out.exists()


# ─────────────────────────────────────────────────────────────────────────────
# to_surfer_xyz
# ─────────────────────────────────────────────────────────────────────────────


def test_to_surfer_xyz_header_and_row_count(tmp_path):
    out = export.to_surfer_xyz(_uniform_model(), tmp_path / "m.dat")
    lines = out.read_text().splitlines()
    assert lines[0] == "X\tY\tZ"
    assert len(lines) == 1 + 3 * 2  # header + n_x * n_z real cells


def test_to_surfer_xyz_no_header(tmp_path):
    out = export.to_surfer_xyz(_uniform_model(), tmp_path / "m.dat", header=False)
    lines = out.read_text().splitlines()
    assert lines[0] != "X\tY\tZ"
    assert len(lines) == 3 * 2


def test_to_surfer_xyz_skips_nan_cells(tmp_path):
    out = export.to_surfer_xyz(_model(), tmp_path / "m.dat")
    text = out.read_text()
    assert "nan" not in text.lower()
    # 2 stations x 3 depths = 6 cells, minus the one real NaN = 5 rows + header
    assert len(text.splitlines()) == 1 + 5


def test_to_surfer_xyz_exact_values_no_resampling(tmp_path):
    out = export.to_surfer_xyz(_uniform_model(), tmp_path / "m.dat")
    rows = out.read_text().splitlines()[1:]
    first = rows[0].split("\t")
    assert first == ["0.000", "-10.000", "2.00000"]


def test_to_surfer_xyz_custom_delimiter(tmp_path):
    out = export.to_surfer_xyz(_uniform_model(), tmp_path / "m.dat", delimiter=",")
    assert "," in out.read_text().splitlines()[0]


def test_to_surfer_xyz_elevation_shifts_y(tmp_path):
    m = _uniform_model()
    flat = export.to_surfer_xyz(m, tmp_path / "flat.dat")
    draped = export.to_surfer_xyz(
        m, tmp_path / "draped.dat", elevation=np.array([50.0, 50.0, 50.0]),
        chainage=m.x_centers,
    )
    y_flat = float(flat.read_text().splitlines()[1].split("\t")[1])
    y_draped = float(draped.read_text().splitlines()[1].split("\t")[1])
    assert y_draped == pytest.approx(y_flat + 50.0)


def test_to_surfer_xyz_accepts_raw_array_triple(tmp_path):
    x = np.array([0.0, 100.0])
    z = np.array([10.0])
    rho = np.array([[2.0, 2.1]])
    out = export.to_surfer_xyz((x, z, rho), tmp_path / "m.dat")
    assert out.exists()


def test_to_surfer_xyz_creates_parent_dirs(tmp_path):
    out = export.to_surfer_xyz(_uniform_model(), tmp_path / "nested" / "m.dat")
    assert out.exists()
