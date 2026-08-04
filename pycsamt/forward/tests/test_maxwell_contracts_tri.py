"""Tests for solver-neutral triangular-mesh Maxwell contracts."""

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    ForwardResult,
    ReceiverSet,
    SolverDiagnostics,
    TriMesh,
    TriProblem,
)


def _mesh():
    # A small quad split into two triangles, plus one extra node/triangle so
    # region_ids/boundary_segments have something non-trivial to describe.
    return TriMesh(
        nodes_m=[[0, 0], [100, 0], [100, 80], [0, 80], [150, 40]],
        triangles=[[0, 1, 2], [0, 2, 3], [1, 4, 2]],
        region_ids=[0, 0, 1],
        boundary_segments=[[0, 1], [1, 4], [4, 2], [2, 3], [3, 0]],
    )


def _problem():
    mesh = _mesh()
    return TriProblem(
        mesh,
        np.full(mesh.shape, 0.01),
        [10, 1],
        ReceiverSet([[50, 0], [120, 0]], ["S00", "S01"]),
        metadata={"model": "half-space"},
    )


def test_mesh_geometry_serialization_and_readonly():
    mesh = _mesh()
    assert mesh.dimension == 2
    assert mesh.n_nodes == 5
    assert mesh.n_triangles == 3
    assert mesh.shape == (3,)
    assert mesh.region_ids.tolist() == [0, 0, 1]
    np.testing.assert_allclose(mesh.triangle_areas_m2, [4000.0, 4000.0, 2000.0])
    assert TriMesh.from_dict(mesh.to_dict()).to_dict() == mesh.to_dict()
    with pytest.raises(ValueError):
        mesh.nodes_m[0, 0] = 2
    with pytest.raises(ValueError):
        mesh.triangles[0, 0] = 4


def test_mesh_rejects_degenerate_and_out_of_range_triangles():
    with pytest.raises(ValueError, match="positive area"):
        TriMesh([[0, 0], [1, 0], [2, 0]], [[0, 1, 2]])
    with pytest.raises(ValueError, match="distinct nodes"):
        TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 0, 1]])
    with pytest.raises(ValueError, match="reference nodes_m"):
        TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 3]])
    with pytest.raises(ValueError, match="one entry per triangle"):
        TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]], region_ids=[0, 1])


def test_mesh_defaults_region_ids_to_single_region():
    mesh = TriMesh([[0, 0], [1, 0], [0, 1]], [[0, 1, 2]])
    assert mesh.region_ids.tolist() == [0]
    assert mesh.boundary_segments is None
    centroid = mesh.triangle_centroids_m
    np.testing.assert_allclose(centroid, [[1 / 3, 1 / 3]])


def test_problem_hash_roundtrip_and_input_sensitivity(tmp_path):
    problem = _problem()
    restored = TriProblem.from_npz(problem.to_npz(tmp_path / "problem.npz"))
    assert restored.problem_hash == problem.problem_hash
    changed = TriProblem(
        problem.mesh,
        problem.conductivity_s_m * 2,
        problem.frequencies_hz,
        problem.receivers,
        problem.components,
    )
    assert changed.problem_hash != problem.problem_hash
    with pytest.raises(ValueError):
        problem.conductivity_s_m[0] = 3


def test_problem_rejects_dimension_components_and_physics_errors():
    mesh = _mesh()
    with pytest.raises(ValueError, match="dimensions"):
        TriProblem(
            mesh,
            np.ones(mesh.shape),
            [1],
            ReceiverSet([[0, 0, 0]], ["S"]),
        )
    with pytest.raises(ValueError, match="only zxy"):
        TriProblem(
            mesh,
            np.ones(mesh.shape),
            [1],
            ReceiverSet([[0, 0]], ["S"]),
            ("zxx",),
        )
    with pytest.raises(ValueError, match="positive"):
        TriProblem(mesh, np.zeros(mesh.shape), [1], ReceiverSet([[0, 0]], ["S"]))
    with pytest.raises(ValueError, match="unique"):
        TriProblem(
            mesh, np.ones(mesh.shape), [1, 1], ReceiverSet([[0, 0]], ["S"])
        )
    with pytest.raises(TypeError, match="TriMesh"):
        TriProblem(
            object(), np.ones(mesh.shape), [1], ReceiverSet([[0, 0]], ["S"])
        )


def test_forward_result_validates_against_tri_problem():
    problem = _problem()
    diagnostics = SolverDiagnostics(
        [[True], [True]], [[3], [4]], [[1e-9], [2e-9]], 0.2
    )
    impedance = np.ones((2, 2, 2), dtype=complex) * (1 + 2j)
    result = ForwardResult(
        problem.problem_hash,
        problem.frequencies_hz,
        problem.receivers.names,
        problem.components,
        impedance,
        None,
        "verified-demo",
        "1.0",
        diagnostics,
    )
    assert result.success
    result.validate_against(problem)
    with pytest.raises(ValueError, match="problem_hash"):
        ForwardResult(
            "0" * 64,
            problem.frequencies_hz,
            problem.receivers.names,
            problem.components,
            impedance,
            None,
            "demo",
            "1",
            diagnostics,
        ).validate_against(problem)


def test_result_npz_roundtrip(tmp_path):
    problem = _problem()
    diagnostics = SolverDiagnostics(
        [[True], [True]], [[0], [0]], [[1e-10], [2e-10]], 0.1
    )
    result = ForwardResult(
        problem.problem_hash,
        problem.frequencies_hz,
        problem.receivers.names,
        problem.components,
        np.full((2, 2, 2), 1 + 1j),
        None,
        "backend",
        "2.0",
        diagnostics,
        {"solver": "direct"},
    )
    restored = ForwardResult.from_npz(result.to_npz(tmp_path / "result.npz"))
    restored.validate_against(problem)
    np.testing.assert_array_equal(restored.impedance_v_a, result.impedance_v_a)
    assert restored.provenance() == result.provenance()
