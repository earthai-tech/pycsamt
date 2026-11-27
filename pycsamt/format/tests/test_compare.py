"""Tests for pycsamt.format.compare (sample_log10, compare_log10)."""

import numpy as np
import pytest

from pycsamt.format import (
    Grid2DGeometry,
    Grid3DGeometry,
    PCSFModel,
    UnstructuredMeshGeometry,
    compare_log10,
    sample_log10,
)


def _grid2d():
    return PCSFModel(
        geometry=Grid2DGeometry(x=[0.0, 10.0, 20.0], z=[5.0, 15.0]),
        resistivity=10.0 ** np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]),
    )


def test_grid2d_linear_in_log10_between_centres():
    out = sample_log10(_grid2d(), [[5.0, 10.0], [15.0, 5.0]])
    np.testing.assert_allclose(out, [1.5, 2.5])


def test_grid2d_edge_half_cell_takes_edge_value_and_outside_is_nan():
    out = sample_log10(_grid2d(), [[0.0, 0.0], [25.0, 10.0], [26.0, 10.0], [10.0, 25.0]])
    np.testing.assert_allclose(out[:2], [1.0, 3.0])
    assert np.isnan(out[2]) and np.isnan(out[3])


def test_grid2d_uses_explicit_nodes_when_present():
    geo = Grid2DGeometry(x=[0.0, 10.0], z=[5.0], x_nodes=[-50.0, 5.0, 15.0],
                         z_nodes=[0.0, 10.0])
    m = PCSFModel(geometry=geo, resistivity=10.0 ** np.array([[1.0, 2.0]]))
    out = sample_log10(m, [[-40.0, 5.0], [-60.0, 5.0]])
    assert out[0] == pytest.approx(1.0)
    assert np.isnan(out[1])


def test_grid3d_trilinear_centre_is_mean():
    m = PCSFModel(
        geometry=Grid3DGeometry(x=[0.0, 10.0], y=[0.0, 10.0], z=[5.0, 15.0]),
        resistivity=10.0 ** np.arange(8.0).reshape(2, 2, 2),
    )
    assert sample_log10(m, [[5.0, 5.0, 10.0]])[0] == pytest.approx(3.5)


def test_mesh_point_location_is_piecewise_constant():
    geo = UnstructuredMeshGeometry(
        nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
        connectivity=np.array([[0, 1, 2], [1, 3, 2]]),
        region_ids=np.array([0, 1]),
    )
    m = PCSFModel(geometry=geo, resistivity=np.array([10.0, 100.0]))
    out = sample_log10(m, [[0.2, 0.2], [0.8, 0.8], [2.0, 2.0]])
    np.testing.assert_allclose(out[:2], [1.0, 2.0])
    assert np.isnan(out[2])


def test_wrong_point_dimension_raises():
    with pytest.raises(ValueError):
        sample_log10(_grid2d(), [[0.0, 0.0, 0.0]])


def test_compare_log10_statistics_ignore_nan():
    res = compare_log10([1.0, 2.0, np.nan, 4.0], [1.5, 2.0, 3.0, 4.5])
    assert res.n == 3
    assert res.bias == pytest.approx(-1.0 / 3.0)
    assert res.rms == pytest.approx(np.sqrt((0.25 + 0.0 + 0.25) / 3.0))
    assert np.isnan(res.delta[2])


def test_compare_log10_shape_mismatch_raises():
    with pytest.raises(ValueError):
        compare_log10([1.0, 2.0], [1.0])


def test_mesh_node_average_geometric_and_recorded():
    from pycsamt.format.adapters import mesh_to_pcsf

    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    conn = np.array([[0, 1, 2]])
    node_rho = np.array([10.0, 100.0, 1000.0])
    arith = mesh_to_pcsf(nodes, conn, resistivity_by_node=node_rho)
    geom = mesh_to_pcsf(nodes, conn, resistivity_by_node=node_rho,
                        node_average="geometric")
    assert arith.resistivity[0] == pytest.approx(370.0)
    assert geom.resistivity[0] == pytest.approx(100.0)
    assert geom.metadata["node_average"] == "geometric"
    with pytest.raises(ValueError):
        mesh_to_pcsf(nodes, conn, resistivity_by_node=node_rho, node_average="median")


def test_log_uncertainty_conversions():
    from pycsamt.format import log_uncertainty

    assert log_uncertainty(0.1, "log10", "ln") == pytest.approx(0.1 * np.log(10))
    assert log_uncertainty(0.1 * np.log(10), "ln", "log10") == pytest.approx(0.1)
    s_lin = log_uncertainty(0.01, "log10", "linear", rho=100.0)
    assert s_lin == pytest.approx(100.0 * np.log(10) * 0.01)
    assert log_uncertainty(s_lin, "linear", "log10", rho=100.0) == pytest.approx(0.01)
    with pytest.raises(ValueError):
        log_uncertainty(0.1, "log10", "linear")


def test_inversion_record_roundtrip_and_check(tmp_path):
    from pycsamt.format import (InversionRecord, attach_inversion_record,
                                check_inversion_record, get_inversion_record,
                                read_pcsf, write_pcsf)

    rec = InversionRecord(
        solver="Occam2D",
        data={"components": ["TE rho"], "frequency_range_hz": [10.0, 1e4],
              "n_stations": 3, "errors": "10 % rho"},
        regularization={"type": "smoothness"},
        misfit={"final_rms": 1.0, "termination": "target_reached"},
    )
    assert check_inversion_record(rec)["required"] == []
    model = attach_inversion_record(_grid2d(), rec)
    path = write_pcsf(model, tmp_path / "m.pcsf")
    back = get_inversion_record(read_pcsf(path))
    assert back == rec
    partial = InversionRecord(solver="ModEM")
    assert "misfit.final_rms" in check_inversion_record(partial)["required"]
    with pytest.raises(ValueError):
        InversionRecord(solver="x", misfit={"termination": "done"})


def test_conductivity_encodings_roundtrip(tmp_path):
    from pycsamt.format import read_pcsf, read_pcsm, write_pcsf, write_pcsm
    from pycsamt.format.adapters import grid3d_to_pcsf

    v = np.array([[[-1.0, -2.0], [np.nan, 0.5]]])
    m = grid3d_to_pcsf(x=[0.0, 1.0], y=[0.0, 1.0], z=[10.0], resistivity=v,
                       encoding="log10_conductivity", crs="EPSG:4326")
    np.testing.assert_allclose(m.resistivity[0, 0], [10.0, 100.0])
    p = write_pcsf(m, tmp_path / "a.pcsf")
    b = read_pcsm(write_pcsm(read_pcsf(p), tmp_path / "a.pcsm"))
    assert b.resistivity_native_encoding == "log10_conductivity"
    assert np.array_equal(b.resistivity_native, v, equal_nan=True)
    s = grid3d_to_pcsf(x=[0.0], y=[0.0], z=[1.0], resistivity=np.array([[[0.01]]]),
                       encoding="conductivity")
    assert s.resistivity[0, 0, 0] == pytest.approx(100.0)
