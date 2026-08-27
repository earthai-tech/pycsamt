# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.format.regrid.mesh_to_grid2d — the regridded
grid2d convenience export for mesh_unstructured models, closing the
gap SPEC.md section 7 / the GMD paper's Discussion previously
documented as deliberately unimplemented.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pycsamt.format.regrid import mesh_to_grid2d
from pycsamt.format.schema import (
    Grid2DGeometry,
    PCSFModel,
    UnstructuredMeshGeometry,
)

triangle = pytest.importorskip("triangle")


def _two_triangle_model() -> PCSFModel:
    # A 1x1 square split into two triangles along the diagonal, real
    # (deterministic) per-triangle resistivity: 10.0 for the lower-left
    # triangle, 20.0 for the upper-right one.
    geometry = UnstructuredMeshGeometry(
        nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
        connectivity=np.array([[0, 1, 2], [1, 3, 2]]),
        region_ids=np.array([0, 1]),
    )
    return PCSFModel(
        geometry=geometry,
        resistivity=np.array([10.0, 20.0]),
        source_backend="mare2dem",
        description="synthetic two-triangle mesh",
    )


class TestMeshToGrid2DSynthetic:
    def test_rejects_non_mesh_model(self):
        model = PCSFModel(
            geometry=Grid2DGeometry(x=np.array([0.0, 1.0]), z=np.array([0.0, 1.0])),
            resistivity=np.ones((2, 2)),
        )
        with pytest.raises(ValueError, match="mesh_unstructured"):
            mesh_to_grid2d(model)

    def test_rejects_region_collapsed_resistivity(self):
        geometry = UnstructuredMeshGeometry(
            nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
            connectivity=np.array([[0, 1, 2], [1, 3, 2]]),
            region_ids=np.array([0, 1]),
        )
        model = PCSFModel(geometry=geometry, resistivity=np.array([10.0]))
        with pytest.raises(ValueError, match="per-cell"):
            mesh_to_grid2d(model)

    def test_rejects_non_xz_plane(self):
        geometry = UnstructuredMeshGeometry(
            nodes=np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]),
            connectivity=np.array([[0, 1, 2], [1, 3, 2]]),
            region_ids=np.array([0, 1]),
            plane="xy",
        )
        model = PCSFModel(geometry=geometry, resistivity=np.array([10.0, 20.0]))
        with pytest.raises(NotImplementedError, match="plane"):
            mesh_to_grid2d(model)

    def test_result_is_grid2d_with_requested_resolution(self):
        model = _two_triangle_model()
        regridded = mesh_to_grid2d(model, nx=9, nz=7)
        assert regridded.kind == "grid2d"
        assert regridded.resistivity.shape == (7, 9)
        regridded.validate()

    def test_interior_points_recover_the_correct_triangle_value(self):
        model = _two_triangle_model()
        regridded = mesh_to_grid2d(model, nx=21, nz=21)
        # Lower-left corner sits deep inside the lower-left triangle
        # (region 10.0); upper-right corner deep inside the other
        # (region 20.0) -- real, deterministic point-location, not an
        # approximation.
        assert regridded.resistivity[0, 0] == pytest.approx(10.0)
        assert regridded.resistivity[-1, -1] == pytest.approx(20.0)

    def test_provenance_is_marked_synthesized(self):
        model = _two_triangle_model()
        regridded = mesh_to_grid2d(model)
        assert regridded.metadata["synthesized"] is True
        assert regridded.metadata["regridded_from"] == "mesh_unstructured"
        assert "synthesized" in regridded.description
        assert regridded.source_backend == "mare2dem"

    def test_no_native_or_region_fields_survive(self):
        model = _two_triangle_model()
        regridded = mesh_to_grid2d(model)
        assert regridded.resistivity_native is None
        assert regridded.resistivity_native_encoding is None
        assert regridded.resistivity_by_region is None

    def test_round_trips_through_pcsf(self, tmp_path):
        from pycsamt.format.io import read_pcsf, write_pcsf

        model = _two_triangle_model()
        regridded = mesh_to_grid2d(model, nx=11, nz=9)
        path = write_pcsf(regridded, tmp_path / "regridded.pcsf")
        restored = read_pcsf(path)
        np.testing.assert_allclose(
            restored.resistivity, regridded.resistivity, equal_nan=True
        )


_DATA_DIR = Path(__file__).parents[3] / "data" / "mare2dem" / "demo_mt_inversion"
_POLY_PATH = _DATA_DIR / "demo.poly"
_SKIP_REAL = pytest.mark.skipif(
    not _POLY_PATH.exists(), reason=f"bundled MARE2DEM data not found: {_POLY_PATH}"
)


def _real_model() -> PCSFModel:
    from pycsamt.forward.maxwell.contracts_tri import TriMesh
    from pycsamt.models.mare2dem.iotools.poly import read_poly
    from pycsamt.models.mare2dem.results import InversionResult
    from pycsamt.format.adapters.mare2dem import mare2dem_to_pcsf

    result = InversionResult(workdir=_DATA_DIR)
    poly = read_poly(_POLY_PATH)
    pslg = {"vertices": poly.nodes, "segments": poly.segments - 1, "regions": poly.regions}
    out = triangle.triangulate(pslg, "pA")
    region_ids = np.round(out["triangle_attributes"].ravel()).astype(np.int64)
    mesh = TriMesh(nodes_m=out["vertices"], triangles=out["triangles"], region_ids=region_ids)
    return mare2dem_to_pcsf(result, mesh, created_by="pytest")


@_SKIP_REAL
class TestMeshToGrid2DRealData:
    def test_regrids_the_real_demo_mt_inversion_mesh(self):
        model = _real_model()
        regridded = mesh_to_grid2d(model, nx=120, nz=80)
        assert regridded.kind == "grid2d"
        assert regridded.resistivity.shape == (80, 120)
        # Most of a real, dense mesh's bounding box is covered -- some
        # nan is expected near the mesh's own irregular boundary, but
        # not everywhere.
        finite = np.isfinite(regridded.resistivity)
        assert finite.mean() > 0.5

    def test_sampled_values_are_real_mesh_resistivities_not_fabricated(self):
        # Every finite value in the regrid must be one of the mesh's
        # own real per-triangle resistivities -- interpolation here is
        # nearest-triangle point-location, not a blend, so no value
        # outside the source set can appear.
        model = _real_model()
        regridded = mesh_to_grid2d(model, nx=60, nz=40)
        finite_values = set(np.round(regridded.resistivity[np.isfinite(regridded.resistivity)], 6))
        real_values = set(np.round(model.resistivity, 6))
        assert finite_values.issubset(real_values)
