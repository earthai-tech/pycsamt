"""Tests for OccamMesh.read()."""

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parents[4] / "data" / "occam2D"
MESH_FILE = DATA_DIR / "Occam2DMesh"

pytestmark = pytest.mark.skipif(
    not MESH_FILE.exists(), reason="sample occam2D data not available"
)


@pytest.fixture(scope="module")
def mesh():
    from pycsamt.models.occam2d.mesh import OccamMesh

    return OccamMesh.read(MESH_FILE)


def test_mesh_n_xcells(mesh):
    assert mesh.n_xcells == 576


def test_mesh_n_zcells(mesh):
    assert mesh.n_zcells == 31


def test_mesh_x_widths_length(mesh):
    assert len(mesh.x_widths) == 576


def test_mesh_z_widths_length(mesh):
    assert len(mesh.z_widths) == 31


def test_mesh_x_widths_first(mesh):
    assert mesh.x_widths[0] == pytest.approx(279.0, rel=1e-4)


def test_mesh_z_widths_first(mesh):
    assert mesh.z_widths[0] == pytest.approx(5.0, rel=1e-4)


def test_mesh_z_widths_last(mesh):
    assert mesh.z_widths[-1] == pytest.approx(3000.0, rel=1e-4)


def test_mesh_x_nodes_length(mesh):
    assert len(mesh.x_nodes) == 577


def test_mesh_z_nodes_length(mesh):
    assert len(mesh.z_nodes) == 32


def test_mesh_x_nodes_starts_at_zero(mesh):
    assert mesh.x_nodes[0] == pytest.approx(0.0)


def test_mesh_z_nodes_starts_at_zero(mesh):
    assert mesh.z_nodes[0] == pytest.approx(0.0)


def test_mesh_x_nodes_monotone(mesh):
    import numpy as np

    assert bool((np.diff(mesh.x_nodes) > 0).all())


def test_mesh_z_nodes_monotone(mesh):
    import numpy as np

    assert bool((np.diff(mesh.z_nodes) > 0).all())


def test_mesh_cell_rows_nonempty(mesh):
    assert len(mesh.cell_rows) > 0


def test_mesh_cell_rows_width(mesh):
    for row in mesh.cell_rows:
        assert len(row) == 576


def test_mesh_all_free_params(mesh):
    for row in mesh.cell_rows:
        assert set(row) == {"?"}


def test_mesh_n_airlayers(mesh):
    assert mesh.n_airlayers == 0


def test_mesh_comment(mesh):
    assert "MESH FILE" in mesh.comment.upper()


def test_mesh_n_params_positive(mesh):
    assert mesh.n_params > 0


# -----------------------------------------------------------------------
# Defaults (no file)
# -----------------------------------------------------------------------


def test_mesh_properties_defaults():
    from pycsamt.models.occam2d.mesh import OccamMesh

    mesh = OccamMesh()
    assert mesh.n_xcells == 0
    assert mesh.n_zcells == 0
    assert mesh.n_params == 0


# -----------------------------------------------------------------------
# Error handling
# -----------------------------------------------------------------------


def test_missing_file_raises():
    from pycsamt.models.occam2d.mesh import OccamMesh

    with pytest.raises(FileNotFoundError):
        OccamMesh.read("/nonexistent/Occam2DMesh")
