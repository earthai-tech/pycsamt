# -*- coding: utf-8 -*-
"""Tests for OccamMesh."""

import pytest
from pathlib import Path

DATA_DIR  = Path(__file__).parents[4] / "data" / "occam2D"
MESH_FILE = DATA_DIR / "Occam2DMesh"


@pytest.mark.skipif(not MESH_FILE.exists(), reason="sample data not available")
def test_read_mesh_file():
    from pycsamt.models.occam2d.mesh import OccamMesh
    mesh = OccamMesh.read(MESH_FILE)
    assert mesh.n_xcells > 0
    assert mesh.n_zcells > 0


def test_mesh_properties_defaults():
    from pycsamt.models.occam2d.mesh import OccamMesh
    mesh = OccamMesh()
    assert mesh.n_xcells == 0
    assert mesh.n_zcells == 0
    assert mesh.n_params == 0
