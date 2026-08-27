# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Phase 7 tests — pcsf_to_point_cloud, the generic 3-D-view extraction
every geometry kind funnels through (no backend-specific glue).

Cross-checked against real data from every earlier phase: Occam2D
(grid2d), ModEM (grid3d), MARE2DEM (mesh_unstructured).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.format import Grid2DGeometry, PCSFModel
from pycsamt.format.multiline import build_multiline_pcsf
from pycsamt.format.pointcloud import pcsf_to_point_cloud

_ROOT = Path(__file__).parents[3]
_OCCAM_DIR = _ROOT / "data" / "occam2D"
_MODEM_DIR = _ROOT / "data" / "modem" / "willy_27freq_watex_line02_sample"
_MARE_DIR = _ROOT / "data" / "mare2dem" / "demo_mt_inversion"


def _grid2d_model() -> PCSFModel:
    geometry = Grid2DGeometry(x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]))
    return PCSFModel(
        geometry=geometry,
        resistivity=np.array([[100.0, 110.0], [50.0, 55.0]]),
    )


class TestGrid2DPointCloud:
    def test_shape_and_z_sign(self):
        cloud = pcsf_to_point_cloud(_grid2d_model())
        assert cloud.x.shape == (4,)
        # z is elevation-like (positive up); depth 10 m -> z == -10.
        assert set(np.round(cloud.z, 3).tolist()) == {-10.0, -50.0}

    def test_value_is_log10_rho(self):
        cloud = pcsf_to_point_cloud(_grid2d_model())
        np.testing.assert_allclose(sorted(cloud.value), sorted(np.log10(
            [100.0, 110.0, 50.0, 55.0]
        )))

    def test_all_y_are_zero_for_a_single_line(self):
        cloud = pcsf_to_point_cloud(_grid2d_model())
        np.testing.assert_array_equal(cloud.y, np.zeros(4))


class TestMultilinePointCloud:
    def test_lines_separated_by_offset_y(self):
        profiles = {
            "L1": {"x": np.array([0.0, 100.0]), "z": np.array([10.0, 50.0]),
                   "rho": np.array([[100.0, 110.0], [50.0, 55.0]])},
            "L2": {"x": np.array([0.0, 100.0]), "z": np.array([10.0, 50.0]),
                   "rho": np.array([[200.0, 210.0], [90.0, 95.0]])},
        }
        model = build_multiline_pcsf(profiles, cache_derived_volume=False)
        cloud = pcsf_to_point_cloud(model)
        assert cloud.x.shape == (8,)
        assert set(np.unique(cloud.y).tolist()) == {0.0, 1000.0}


class TestSubsampling:
    def test_caps_at_max_points_reproducibly(self):
        geometry = Grid2DGeometry(
            x=np.linspace(0, 1000, 200), z=np.linspace(1, 500, 200)
        )
        model = PCSFModel(
            geometry=geometry, resistivity=np.full((200, 200), 100.0)
        )
        cloud1 = pcsf_to_point_cloud(model, max_points=1000, seed=42)
        cloud2 = pcsf_to_point_cloud(model, max_points=1000, seed=42)
        assert cloud1.x.size == 1000
        assert "subsampled" in cloud1.label
        np.testing.assert_array_equal(cloud1.x, cloud2.x)

    def test_no_subsampling_when_under_the_cap(self):
        cloud = pcsf_to_point_cloud(_grid2d_model(), max_points=1000)
        assert cloud.x.size == 4
        assert "subsampled" not in cloud.label


class TestNonFiniteDropped:
    def test_nan_cells_are_dropped_not_zeroed(self):
        geometry = Grid2DGeometry(x=np.array([0.0, 100.0]), z=np.array([10.0, 50.0]))
        model = PCSFModel(
            geometry=geometry,
            resistivity=np.array([[100.0, np.nan], [50.0, 55.0]]),
        )
        cloud = pcsf_to_point_cloud(model)
        assert cloud.x.shape == (3,)


class TestErrors:
    def test_rejects_mesh_unstructured_with_region_shaped_resistivity(self):
        from pycsamt.format.schema import UnstructuredMeshGeometry

        # 3 triangles sharing only 2 distinct region ids, so the
        # region-count (2) and triangle-count (3) genuinely differ --
        # region-shaped resistivity must not silently line up by luck.
        geometry = UnstructuredMeshGeometry(
            nodes=np.array(
                [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.5, 1.5]]
            ),
            connectivity=np.array([[0, 1, 2], [1, 3, 2], [2, 3, 4]]),
            region_ids=np.array([1, 2, 1]),
        )
        model = PCSFModel(geometry=geometry, resistivity=np.array([10.0, 20.0]))
        with pytest.raises(ValueError, match="per-cell"):
            pcsf_to_point_cloud(model)


# ---------------------------------------------------------------------
# Real data, every backend already validated in earlier phases
# ---------------------------------------------------------------------


@pytest.mark.skipif(not _OCCAM_DIR.exists(), reason="bundled occam2D data missing")
def test_real_occam2d_grid2d():
    from pycsamt.format.adapters.occam2d import occam2d_to_pcsf
    from pycsamt.models.occam2d.results import InversionResult

    model = occam2d_to_pcsf(InversionResult(workdir=_OCCAM_DIR))
    cloud = pcsf_to_point_cloud(model)
    assert cloud.x.size > 0
    assert np.all(np.isfinite(cloud.value))
    assert cloud.z.max() <= 0.0  # all at/below the surface


@pytest.mark.skipif(not _MODEM_DIR.exists(), reason="bundled ModEM data missing")
def test_real_modem_grid3d_respects_max_points():
    from pycsamt.format.adapters.modem3d import modem3d_to_pcsf
    from pycsamt.models.modem.results import InversionResult

    result = InversionResult(workdir=_MODEM_DIR, load_data=True)
    model = modem3d_to_pcsf(result)
    cloud = pcsf_to_point_cloud(model, max_points=20_000)
    assert cloud.x.size == 20_000
    assert "grid3d" in cloud.label


@pytest.mark.skipif(not _MARE_DIR.exists(), reason="bundled MARE2DEM data missing")
def test_real_mare2dem_mesh_unstructured():
    triangle = pytest.importorskip("triangle")
    from pycsamt.forward.maxwell.contracts_tri import TriMesh
    from pycsamt.format.adapters.mare2dem import mare2dem_to_pcsf
    from pycsamt.models.mare2dem.iotools.poly import read_poly
    from pycsamt.models.mare2dem.results import InversionResult

    pf = read_poly(_MARE_DIR / "demo.poly")
    pslg = {"vertices": pf.nodes, "segments": pf.segments - 1, "regions": pf.regions}
    out = triangle.triangulate(pslg, "pA")
    region_ids = np.round(out["triangle_attributes"].ravel()).astype(np.int64)
    mesh = TriMesh(
        nodes_m=out["vertices"], triangles=out["triangles"], region_ids=region_ids
    )
    model = mare2dem_to_pcsf(InversionResult(workdir=_MARE_DIR), mesh)
    cloud = pcsf_to_point_cloud(model, max_points=5000)
    assert cloud.x.size == 5000
    assert np.all(np.isfinite(cloud.value))
