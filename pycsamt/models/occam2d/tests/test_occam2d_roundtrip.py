# -*- coding: utf-8 -*-
"""Round-trip tests: read → write → read-back → compare."""

import math
import pytest
from pathlib import Path

import numpy as np

DATA_DIR    = Path(__file__).parents[4] / "data" / "occam2D"
STARTUP_FILE = DATA_DIR / "Startup"
ITER_FILE    = DATA_DIR / "ITER17.iter"
MESH_FILE    = DATA_DIR / "Occam2DMesh"
MODEL_FILE   = DATA_DIR / "Occam2DModel"
DATA_FILE    = DATA_DIR / "OccamDataFile.dat"

pytestmark = pytest.mark.skipif(
    not DATA_DIR.exists(), reason="sample occam2D data not available"
)


# -----------------------------------------------------------------------
# OccamStartup round-trip
# -----------------------------------------------------------------------

def test_startup_roundtrip(tmp_path):
    from pycsamt.models.occam2d.startup import OccamStartup
    orig = OccamStartup.read(STARTUP_FILE)
    out  = tmp_path / "Startup_rt"
    orig.write(out)
    back = OccamStartup.read(out)

    assert back.format_str      == orig.format_str
    assert back.model_file      == orig.model_file
    assert back.data_file       == orig.data_file
    assert back.max_iterations  == orig.max_iterations
    assert back.roughness_type  == orig.roughness_type
    assert back.iteration       == 0
    assert back.n_params        == orig.n_params
    assert math.isclose(back.target_misfit,  orig.target_misfit,  rel_tol=1e-5)
    assert math.isclose(back.lagrange_value, orig.lagrange_value, rel_tol=1e-4)
    assert np.allclose(back.param_values, orig.param_values, atol=1e-4)


# -----------------------------------------------------------------------
# OccamIter round-trip
# -----------------------------------------------------------------------

def test_iter_roundtrip(tmp_path):
    from pycsamt.models.occam2d.startup import OccamIter
    orig = OccamIter.read(ITER_FILE)
    out  = tmp_path / "ITER17_rt.iter"
    orig.write(out)
    back = OccamIter.read(out)

    assert back.iteration  == orig.iteration
    assert back.n_params   == orig.n_params
    assert math.isclose(back.misfit_value,   orig.misfit_value,   rel_tol=1e-4)
    assert math.isclose(back.lagrange_value, orig.lagrange_value, rel_tol=1e-4)
    assert math.isclose(back.roughness_value, orig.roughness_value, rel_tol=1e-4)
    assert back.misfit_reached == orig.misfit_reached
    assert np.allclose(back.param_values, orig.param_values, atol=1e-4)


# -----------------------------------------------------------------------
# OccamMesh round-trip
# -----------------------------------------------------------------------

def test_mesh_roundtrip(tmp_path):
    from pycsamt.models.occam2d.mesh import OccamMesh
    orig = OccamMesh.read(MESH_FILE)
    out  = tmp_path / "Occam2DMesh_rt"
    orig.write(out)
    back = OccamMesh.read(out)

    assert back.n_xcells    == orig.n_xcells
    assert back.n_zcells    == orig.n_zcells
    assert back.n_airlayers == orig.n_airlayers
    assert np.allclose(back.x_widths, orig.x_widths, rtol=1e-4)
    assert np.allclose(back.z_widths, orig.z_widths, rtol=1e-4)
    assert len(back.cell_rows) == len(orig.cell_rows)
    assert back.cell_rows[0] == orig.cell_rows[0]


# -----------------------------------------------------------------------
# OccamModel round-trip
# -----------------------------------------------------------------------

def test_model_roundtrip(tmp_path):
    from pycsamt.models.occam2d.model import OccamModel
    orig = OccamModel.read(MODEL_FILE)
    out  = tmp_path / "Occam2DModel_rt"
    orig.write(out)
    back = OccamModel.read(out)

    assert back.n_layers      == orig.n_layers
    assert back.n_params      == orig.n_params
    assert back.n_exceptions  == orig.n_exceptions
    assert math.isclose(back.binding_offset, orig.binding_offset, abs_tol=1e-6)
    assert back.mesh_file     == orig.mesh_file
    assert back.mesh_type     == orig.mesh_type
    # Verify every layer's params array
    for i, (bl, ol) in enumerate(zip(back.layers, orig.layers)):
        assert bl["n_merge"] == ol["n_merge"], f"layer {i} n_merge mismatch"
        assert bl["n_cols"]  == ol["n_cols"],  f"layer {i} n_cols mismatch"
        assert np.array_equal(bl["params"], ol["params"]), f"layer {i} params mismatch"


# -----------------------------------------------------------------------
# OccamData round-trip
# -----------------------------------------------------------------------

def test_data_roundtrip(tmp_path):
    from pycsamt.models.occam2d.data import OccamData
    orig = OccamData.read(DATA_FILE)
    out  = tmp_path / "OccamDataFile_rt.dat"
    orig.write(out)
    back = OccamData.read(out)

    assert back.n_sites       == orig.n_sites
    assert back.n_frequencies == orig.n_frequencies
    assert back.n_data        == orig.n_data
    assert back.sites         == orig.sites
    assert np.allclose(back.offsets,     orig.offsets,     rtol=1e-4)
    assert np.allclose(back.frequencies, orig.frequencies, rtol=1e-6)
    assert np.allclose(back.data_blocks, orig.data_blocks, rtol=1e-3)
