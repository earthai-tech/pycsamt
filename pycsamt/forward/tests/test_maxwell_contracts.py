"""Tests for solver-neutral Maxwell contracts."""

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    ForwardResult,
    MaxwellMesh,
    MaxwellProblem,
    ReceiverSet,
    SolverDiagnostics,
)


def _problem(dimension=2):
    mesh = (
        MaxwellMesh([0, 100, 250], [0, 50, 150])
        if dimension == 2
        else MaxwellMesh([0, 100, 250], [0, 50, 150], [0, 80, 200])
    )
    coordinates = [[50, 0]] if dimension == 2 else [[50, 40, 0]]
    components = ("zxy", "zyx") if dimension == 2 else ("zxx", "zxy", "zyx", "zyy")
    return MaxwellProblem(
        mesh,
        np.full(mesh.shape, 0.01),
        [10, 1],
        ReceiverSet(coordinates, ["S00"]),
        components,
        metadata={"model": "half-space"},
    )


def test_mesh_geometry_serialization_and_readonly():
    mesh = _problem().mesh
    assert mesh.shape == (2, 2)
    assert mesh.cell_widths_m["x"].tolist() == [100, 150]
    assert MaxwellMesh.from_dict(mesh.to_dict()).to_dict() == mesh.to_dict()
    with pytest.raises(ValueError):
        mesh.x_edges_m[0] = 2
    with pytest.raises(ValueError, match="strictly increasing"):
        MaxwellMesh([0, 1, 1], [0, 1, 2])


def test_receivers_validate_names_dimension_and_orientation():
    receivers = ReceiverSet([[0, 0], [1, 0]], ["A", "B"], 370)
    assert receivers.orientation_deg == 10
    assert ReceiverSet.from_dict(receivers.to_dict()).names == receivers.names
    with pytest.raises(ValueError, match="unique"):
        ReceiverSet([[0, 0], [1, 0]], ["A", "A"])


def test_problem_hash_roundtrip_and_input_sensitivity(tmp_path):
    problem = _problem(3)
    restored = MaxwellProblem.from_npz(problem.to_npz(tmp_path / "problem.npz"))
    assert restored.problem_hash == problem.problem_hash
    changed = MaxwellProblem(
        problem.mesh,
        problem.conductivity_s_m * 2,
        problem.frequencies_hz,
        problem.receivers,
        problem.components,
    )
    assert changed.problem_hash != problem.problem_hash
    with pytest.raises(ValueError):
        problem.conductivity_s_m[0, 0, 0] = 3


def test_problem_rejects_dimension_components_and_physics_errors():
    mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
    with pytest.raises(ValueError, match="dimensions"):
        MaxwellProblem(mesh, np.ones(mesh.shape), [1], ReceiverSet([[0, 0, 0]], ["S"]))
    with pytest.raises(ValueError, match="only zxy"):
        MaxwellProblem(
            mesh,
            np.ones(mesh.shape),
            [1],
            ReceiverSet([[0, 0]], ["S"]),
            ("zxx",),
        )
    with pytest.raises(ValueError, match="positive"):
        MaxwellProblem(mesh, np.zeros(mesh.shape), [1], ReceiverSet([[0, 0]], ["S"]))
    with pytest.raises(ValueError, match="unique"):
        MaxwellProblem(mesh, np.ones(mesh.shape), [1, 1], ReceiverSet([[0, 0]], ["S"]))


def test_diagnostics_and_result_validation():
    problem = _problem()
    diagnostics = SolverDiagnostics([[True], [True]], [[3], [4]], [[1e-9], [2e-9]], 0.2)
    impedance = np.ones((1, 2, 2), dtype=complex) * (1 + 2j)
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
        np.full((1, 2, 2), 1 + 1j),
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


def test_result_masks_nonfinite_values_and_rejects_false_validity():
    diagnostics = SolverDiagnostics([[False]], [[10]], [[0.1]], 1, ("failed",))
    result = ForwardResult(
        "0" * 64,
        [1],
        ["S"],
        ["zxy"],
        [[[complex(np.nan, 0)]]],
        None,
        "demo",
        "1",
        diagnostics,
    )
    assert not result.valid.item()
    assert not result.success
    with pytest.raises(ValueError, match="must be finite"):
        ForwardResult(
            "0" * 64,
            [1],
            ["S"],
            ["zxy"],
            [[[complex(np.nan, 0)]]],
            [[[True]]],
            "demo",
            "1",
            diagnostics,
        )


def test_invalid_diagnostics_are_rejected():
    with pytest.raises(ValueError, match="identical shapes"):
        SolverDiagnostics([[True]], [[1, 2]], [[0.1]], 0)
    with pytest.raises(ValueError, match="non-negative integers"):
        SolverDiagnostics([[True]], [[-1]], [[0.1]], 0)
