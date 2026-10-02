# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Phase 4 tests — the MARE2DEM adapter, against the real bundled
``data/mare2dem/demo_mt_inversion`` dataset (6540 real regions).

The mesh paired with the real ``.resistivity`` table is built with a
genuine, in-process constrained triangulation (the ``triangle`` Python
package — already a hard pycsamt dependency, no external Triangle
binary needed) of the run's own real ``demo.poly`` PSLG, not a
fabricated stand-in: it reproduces the real region partition exactly
(6540/6540 unique region ids, no gaps), confirming ``demo.poly`` really
is the PSLG MARE2DEM triangulated for this run.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

triangle = pytest.importorskip("triangle")

from pycsamt.format import read_pcsf, write_pcsf
from pycsamt.format.adapters.mare2dem import (
    _expand_region_resistivity,
    mare2dem_to_pcsf,
)
from pycsamt.forward.maxwell.contracts_tri import TriMesh
from pycsamt.models.mare2dem.iotools.poly import read_poly
from pycsamt.models.mare2dem.results import InversionResult

_DATA_DIR = Path(__file__).parents[3] / "data" / "mare2dem" / "demo_mt_inversion"
_POLY_PATH = _DATA_DIR / "demo.poly"
_SKIP = pytest.mark.skipif(
    not _POLY_PATH.exists(), reason=f"bundled MARE2DEM data not found: {_POLY_PATH}"
)


def _real_mesh() -> TriMesh:
    pf = read_poly(_POLY_PATH)
    pslg = {
        "vertices": pf.nodes,
        "segments": pf.segments - 1,
        "regions": pf.regions,
    }
    out = triangle.triangulate(pslg, "pA")
    region_ids = np.round(out["triangle_attributes"].ravel()).astype(np.int64)
    return TriMesh(
        nodes_m=out["vertices"], triangles=out["triangles"], region_ids=region_ids
    )


@pytest.fixture(scope="module")
def result() -> InversionResult:
    return InversionResult(workdir=_DATA_DIR)


@pytest.fixture(scope="module")
def mesh() -> TriMesh:
    return _real_mesh()


@pytest.fixture(scope="module")
def model(result, mesh):
    return mare2dem_to_pcsf(
        result,
        mesh,
        created_by="pytest",
        description="MARE2DEM demo_mt_inversion",
    )


@_SKIP
class TestMare2DEMAdapter:
    def test_real_mesh_covers_every_region_exactly(self, result, mesh):
        # Confirms demo.poly really is the PSLG this run's own
        # .resistivity table partitions — not an approximation.
        ids = set(np.unique(mesh.region_ids).tolist())
        assert ids == set(range(1, result.model.num_regions + 1))

    def test_geometry_kind_and_shapes(self, model, mesh):
        assert model.kind == "mesh_unstructured"
        assert model.geometry.connectivity.shape == mesh.triangles.shape
        assert model.resistivity.shape == (mesh.n_triangles,)
        model.validate()

    def test_resistivity_by_region_matches_source_file(self, result, model):
        np.testing.assert_array_equal(
            model.resistivity_by_region, result.model.resistivity[:, 0]
        )

    def test_expanded_resistivity_matches_region_lookup(self, result, model, mesh):
        expected = result.model.resistivity[:, 0][mesh.region_ids - 1]
        np.testing.assert_array_equal(model.resistivity, expected)

    def test_resistivity_is_already_linear_no_native_conversion(self, model):
        # Unlike Occam2D (log10) / ModEM (ln), MARE2DEM's own file is
        # already linear ohm.m -- nothing to preserve as "native".
        assert model.resistivity_native is None
        assert model.resistivity_native_encoding is None

    def test_metadata_carries_provenance(self, result, model):
        assert model.metadata["n_regions"] == result.model.num_regions
        assert model.metadata["anisotropy"] == "isotropic"
        assert model.metadata["converged"] == result.converged

    def test_full_round_trip_through_pcsf_file_is_bit_exact(self, model, tmp_path):
        path = write_pcsf(model, tmp_path / "mare2dem_real.pcsf")
        restored = read_pcsf(path)

        assert restored.source_backend == "mare2dem"
        np.testing.assert_array_equal(restored.resistivity, model.resistivity)
        np.testing.assert_array_equal(
            restored.resistivity_by_region, model.resistivity_by_region
        )
        np.testing.assert_array_equal(
            restored.geometry.connectivity, model.geometry.connectivity
        )
        np.testing.assert_array_equal(
            restored.geometry.region_ids, model.geometry.region_ids
        )


def test_expand_region_resistivity_rejects_out_of_range_id():
    with pytest.raises(ValueError, match="outside"):
        _expand_region_resistivity(
            np.array([10.0, 20.0]), np.array([1, 2, 3])
        )


def test_expand_region_resistivity_basic_lookup():
    result = _expand_region_resistivity(
        np.array([10.0, 20.0, 30.0]), np.array([1, 3, 2, 1])
    )
    np.testing.assert_array_equal(result, [10.0, 30.0, 20.0, 10.0])


def test_raises_when_result_has_no_model():
    empty = InversionResult.__new__(InversionResult)
    empty.model = None
    mesh = TriMesh(
        nodes_m=[[0, 0], [1, 0], [0, 1]],
        triangles=[[0, 1, 2]],
        region_ids=np.array([1]),
    )
    with pytest.raises(ValueError, match="resistivity model"):
        mare2dem_to_pcsf(empty, mesh)


def test_raises_when_mesh_has_no_region_ids():
    class _FakeRF:
        resistivity = np.array([[10.0]])
        anisotropy = "isotropic"

    empty = InversionResult.__new__(InversionResult)
    empty.model = _FakeRF()
    mesh_no_regions = TriMesh(
        nodes_m=[[0, 0], [1, 0], [0, 1]], triangles=[[0, 1, 2]]
    )
    with pytest.raises(ValueError, match="region_ids"):
        mare2dem_to_pcsf(empty, mesh_no_regions)


def test_raises_for_anisotropic_model():
    class _FakeRF:
        resistivity = np.array([[10.0, 1.0]])
        anisotropy = "triaxial"

    empty = InversionResult.__new__(InversionResult)
    empty.model = _FakeRF()
    mesh = TriMesh(
        nodes_m=[[0, 0], [1, 0], [0, 1]],
        triangles=[[0, 1, 2]],
        region_ids=np.array([1]),
    )
    with pytest.raises(ValueError, match="isotropic"):
        mare2dem_to_pcsf(empty, mesh)


# ---------------------------------------------------------------------------
# Fully self-contained tests (no bundled data/mare2dem/ dependency, so the
# coverage gain is real in CI — see test class above which entirely skips
# there since data/mare2dem/ is gitignored).
# ---------------------------------------------------------------------------


class _FakeRF:
    resistivity = np.array([[10.0], [20.0], [30.0]])
    anisotropy = "isotropic"
    num_regions = 3


class _FakeIterRecord:
    def __init__(self, iteration, rms, roughness, lambda_):
        self.iteration = iteration
        self.rms = rms
        self.roughness = roughness
        self.lambda_ = lambda_


class _FakeLog:
    def __init__(self, iterations):
        self.iterations = iterations
        self.final_rms = iterations[-1].rms if iterations else None
        self.n_iterations = len(iterations)
        self.converged = True


def _fake_result(with_log=True):
    result = InversionResult.__new__(InversionResult)
    result.model = _FakeRF()
    result.workdir = Path("fake_workdir")
    result.log = (
        _FakeLog(
            [
                _FakeIterRecord(0, 5.0, 2.0, 100.0),
                _FakeIterRecord(1, 1.2, 1.5, 10.0),
            ]
        )
        if with_log
        else None
    )
    return result


def _fake_mesh():
    return TriMesh(
        nodes_m=[[0, 0], [1, 0], [0, 1], [1, 1]],
        triangles=[[0, 1, 2], [1, 3, 2]],
        region_ids=np.array([1, 3]),
    )


def test_mare2dem_to_pcsf_full_synthetic_build():
    result = _fake_result()
    mesh = _fake_mesh()
    model = mare2dem_to_pcsf(
        result,
        mesh,
        created_by="pytest",
        crs="EPSG:32650",
        description="synthetic MARE2DEM run",
    )

    assert model.kind == "mesh_unstructured"
    assert model.source_backend == "mare2dem"
    assert model.created_by == "pytest"
    assert model.crs == "EPSG:32650"
    assert model.description == "synthetic MARE2DEM run"
    np.testing.assert_array_equal(model.resistivity, [10.0, 30.0])
    np.testing.assert_array_equal(
        model.resistivity_by_region, [10.0, 20.0, 30.0]
    )
    assert model.geometry.plane == "xz"
    np.testing.assert_array_equal(
        model.geometry.connectivity, mesh.triangles
    )

    assert model.metadata["workdir"] == str(Path("fake_workdir"))
    assert model.metadata["final_rms"] == 1.2
    assert model.metadata["n_iterations"] == 2
    assert model.metadata["converged"] is True
    assert model.metadata["n_regions"] == 3
    assert model.metadata["anisotropy"] == "isotropic"

    assert set(model.history) == {"iteration", "rms", "roughness", "lambda"}
    np.testing.assert_array_equal(model.history["iteration"], [0.0, 1.0])
    np.testing.assert_array_equal(model.history["rms"], [5.0, 1.2])
    np.testing.assert_array_equal(model.history["roughness"], [2.0, 1.5])
    np.testing.assert_array_equal(model.history["lambda"], [100.0, 10.0])


def test_mare2dem_to_pcsf_without_log_history_is_empty():
    result = _fake_result(with_log=False)
    mesh = _fake_mesh()
    model = mare2dem_to_pcsf(result, mesh)
    assert model.history == {}


def test_mare2dem_to_pcsf_with_stations_and_dict_survey():
    from pycsamt.format.schema import StationTable

    result = _fake_result()
    mesh = _fake_mesh()
    stations = StationTable(
        name=["S1", "S2"],
        x=np.array([0.0, 10.0]),
        y=np.array([0.0, 0.0]),
        z=np.array([0.0, 0.0]),
    )
    model = mare2dem_to_pcsf(
        result, mesh, stations=stations, survey={"name": "demo survey"}
    )
    assert model.stations is stations
    assert model.survey == {"name": "demo survey"}


def test_mare2dem_to_pcsf_with_survey_object_using_to_dict():
    class _Survey:
        def to_dict(self):
            return {"name": "obj survey", "crew": "A"}

    result = _fake_result()
    mesh = _fake_mesh()
    model = mare2dem_to_pcsf(result, mesh, survey=_Survey())
    assert model.survey == {"name": "obj survey", "crew": "A"}


def test_survey_to_dict_none_and_mapping():
    from pycsamt.format.adapters.mare2dem import _survey_to_dict

    assert _survey_to_dict(None) == {}
    assert _survey_to_dict({"a": 1}) == {"a": 1}
