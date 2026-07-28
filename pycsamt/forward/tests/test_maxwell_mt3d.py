"""Tests for the research-only 3-D MT finite-difference adapter.

Two layers of validation are used, matching the module's own claims:

1. The discrete curl operators are checked independently of any
   physics (uniform field -> zero, linear field -> exact analytic
   curl, and the topological identity ``div(curl(E)) == 0`` as an
   exact matrix product) -- this is what actually caught the sign
   convention during development, before any physics benchmark ran.
2. The assembled solver is checked against the project's own analytic
   benchmarks (half-space passes; layered-earth does not, and that
   failure is asserted explicitly rather than ignored).
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    IncompatibleProblemError,
    MaxwellMesh,
    MaxwellProblem,
    ReceiverSet,
    half_space_benchmark,
    layered_earth_benchmark,
)
from pycsamt.forward.maxwell.mt3d import (
    _MU0,
    MT3DAdapter,
    _curl_e2h,
    _curl_h2e,
    _edge_shapes,
)


def _calibrated_mesh():
    nx, ny, nz = 6, 6, 10
    hx, hy, hz = 200.0, 200.0, 150.0
    mesh = MaxwellMesh(
        np.arange(nx + 1) * hx,
        np.arange(nz + 1) * hz,
        np.arange(ny + 1) * hy,
    )
    receivers = ReceiverSet([[nx * hx / 2, ny * hy / 2, 0.0]], ["S00"])
    return mesh, receivers


def _problem(**changes):
    mesh = MaxwellMesh([0, 1000, 2000], [0, 500, 1000], [0, 1000, 2000])
    values = dict(
        mesh=mesh,
        conductivity_s_m=np.full(mesh.shape, 0.01),
        frequencies_hz=[1.0],
        receivers=ReceiverSet([[1000.0, 1000.0, 0.0]], ["S00"]),
        components=("zxx", "zxy", "zyx", "zyy"),
    )
    values.update(changes)
    return MaxwellProblem(**values)


# --------------------------------------------------------------------------- #
# Discrete curl operators: pure linear-algebra checks, no physics involved
# --------------------------------------------------------------------------- #


def test_curl_e2h_of_uniform_field_is_exactly_zero():
    nx, ny, nz = 3, 4, 5
    hx, hy, hz = 10.0, 20.0, 30.0
    s = _edge_shapes(nx, ny, nz)
    c1 = _curl_e2h(nx, ny, nz, hx, hy, hz)
    e = np.concatenate(
        [np.full(np.prod(s[name]), 1 + 2j) for name in ("ex", "ey", "ez")]
    )
    assert np.max(np.abs(c1 @ e)) == 0.0


def test_curl_e2h_matches_known_analytic_curl_of_a_linear_field():
    # E = (0, 0, c*x): curl(E) = (0, -c, 0) exactly.
    nx, ny, nz = 3, 4, 5
    hx, hy, hz = 10.0, 20.0, 30.0
    s = _edge_shapes(nx, ny, nz)
    c1 = _curl_e2h(nx, ny, nz, hx, hy, hz)
    c_val = 3.0
    x_edges = np.arange(nx + 1) * hx
    ez = np.zeros(s["ez"])
    for i in range(nx + 1):
        ez[i, :, :] = c_val * x_edges[i]
    e = np.concatenate([np.zeros(np.prod(s["ex"])), np.zeros(np.prod(s["ey"])), ez.ravel()])
    h = c1 @ e
    n_hx, n_hy = np.prod(s["hx"]), np.prod(s["hy"])
    hx_arr = h[:n_hx].reshape(s["hx"])
    hy_arr = h[n_hx : n_hx + n_hy].reshape(s["hy"])
    hz_arr = h[n_hx + n_hy :].reshape(s["hz"])
    assert np.max(np.abs(hx_arr)) == 0.0
    assert np.allclose(hy_arr, -c_val)
    assert np.max(np.abs(hz_arr)) == 0.0


def test_curl_h2e_matches_known_analytic_curl_at_interior_edges():
    # H = (0, 0, c*x_center): curl(H)_y = -c at interior Ey edges.
    nx, ny, nz = 3, 4, 5
    hx, hy, hz = 10.0, 20.0, 30.0
    s = _edge_shapes(nx, ny, nz)
    c2 = _curl_h2e(nx, ny, nz, hx, hy, hz)
    c_val = 5.0
    x_centers = (np.arange(nx) + 0.5) * hx
    hz_arr = np.zeros(s["hz"])
    for i in range(nx):
        hz_arr[i, :, :] = c_val * x_centers[i]
    h = np.concatenate(
        [np.zeros(np.prod(s["hx"])), np.zeros(np.prod(s["hy"])), hz_arr.ravel()]
    )
    e = c2 @ h
    n_ex, n_ey = np.prod(s["ex"]), np.prod(s["ey"])
    ey_arr = e[n_ex : n_ex + n_ey].reshape(s["ey"])
    interior = ey_arr[1:nx, :, 1:nz]
    assert np.allclose(interior, -c_val)


def test_div_curl_e_identity_holds_exactly_as_a_matrix():
    """div(curl(E)) == 0 is a topological identity for ANY E, not just
    hand-picked test fields; checking it as an exact matrix product
    catches indexing bugs that special-case fields could miss."""
    nx, ny, nz = 3, 4, 5
    hx, hy, hz = 10.0, 20.0, 30.0
    s = _edge_shapes(nx, ny, nz)
    c1 = _curl_e2h(nx, ny, nz, hx, hy, hz)
    n_hx, n_hy, n_hz = np.prod(s["hx"]), np.prod(s["hy"]), np.prod(s["hz"])

    rows, cols, vals = [], [], []

    def add(r, c, v):
        rows.append(r)
        cols.append(c)
        vals.append(v)

    off_hy, off_hz = n_hx, n_hx + n_hy

    def hx_idx(i, j, k):
        return np.ravel_multi_index((i, j, k), s["hx"])

    def hy_idx(i, j, k):
        return np.ravel_multi_index((i, j, k), s["hy"])

    def hz_idx(i, j, k):
        return np.ravel_multi_index((i, j, k), s["hz"])

    n_cells = nx * ny * nz
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                r = np.ravel_multi_index((i, j, k), (nx, ny, nz))
                add(r, hx_idx(i + 1, j, k), 1.0 / hx)
                add(r, hx_idx(i, j, k), -1.0 / hx)
                add(r, off_hy + hy_idx(i, j + 1, k), 1.0 / hy)
                add(r, off_hy + hy_idx(i, j, k), -1.0 / hy)
                add(r, off_hz + hz_idx(i, j, k + 1), 1.0 / hz)
                add(r, off_hz + hz_idx(i, j, k), -1.0 / hz)

    from scipy import sparse

    n_h = n_hx + n_hy + n_hz
    div = sparse.csr_matrix((vals, (rows, cols)), shape=(n_cells, n_h))
    assert np.max(np.abs((div @ c1).toarray())) == 0.0


# --------------------------------------------------------------------------- #
# Capability declaration and assess() preflight
# --------------------------------------------------------------------------- #


def test_capabilities_declare_3d_full_tensor_and_a_cell_ceiling():
    adapter = MT3DAdapter()
    cap = adapter.capabilities
    assert cap.name == "mt3d"
    assert cap.dimensions == (3,)
    assert set(cap.components) == {"zxx", "zxy", "zyx", "zyy"}
    assert cap.supports_nonuniform_mesh is False
    assert cap.maximum_cells is not None


def test_assess_rejects_buried_receivers():
    problem = _problem(
        receivers=ReceiverSet([[1000.0, 1000.0, 50.0]], ["S00"])
    )
    report = MT3DAdapter().assess(problem)
    assert not report.compatible
    assert any("surface" in message for message in report.errors)


def test_assess_rejects_receivers_outside_horizontal_bounds():
    problem = _problem(
        receivers=ReceiverSet([[50_000.0, 1000.0, 0.0]], ["S00"])
    )
    report = MT3DAdapter().assess(problem)
    assert not report.compatible
    assert any("mesh" in message for message in report.errors)


def test_assess_rejects_non_vacuum_permeability():
    problem = _problem(magnetic_permeability_h_m=_MU0 * 2.0)
    report = MT3DAdapter().assess(problem)
    assert not report.compatible
    assert any("permeability" in message for message in report.errors)


def test_assess_rejects_nonuniform_mesh():
    mesh = MaxwellMesh([0, 1000, 3000], [0, 500, 1000], [0, 1000, 2000])
    problem = _problem(
        mesh=mesh, conductivity_s_m=np.full(mesh.shape, 0.01)
    )
    with pytest.raises(IncompatibleProblemError, match="nonuniform"):
        MT3DAdapter().solve(problem)


def test_assess_rejects_problems_above_the_cell_ceiling():
    problem = _problem()
    small = MT3DAdapter(max_cells=1)
    with pytest.raises(IncompatibleProblemError, match="cell count"):
        small.solve(problem)


def test_incompatible_2d_problem_is_rejected_before_solving():
    mesh = MaxwellMesh([0, 1, 2], [0, 1, 2])
    problem = MaxwellProblem(
        mesh,
        np.full(mesh.shape, 0.01),
        [1.0],
        ReceiverSet([[0.5, 0.0]], ["S"]),
        ("zxy",),
    )
    with pytest.raises(IncompatibleProblemError, match="2-D"):
        MT3DAdapter().solve(problem)


def test_inactive_cells_are_rejected():
    mesh = MaxwellMesh([0, 1000, 2000], [0, 500, 1000], [0, 1000, 2000])
    active = np.ones(mesh.shape, dtype=bool)
    active[0, 0, 0] = False
    problem = _problem(mesh=mesh, active_cells=active)
    with pytest.raises(IncompatibleProblemError, match="inactive"):
        MT3DAdapter().solve(problem)


# --------------------------------------------------------------------------- #
# Solve round-trip
# --------------------------------------------------------------------------- #


def test_solve_returns_canonical_shape_and_succeeds():
    problem = _problem()
    result = MT3DAdapter().solve(problem)
    assert result.shape == (1, 1, 4)
    assert result.success
    assert result.backend_name == "mt3d"
    result.validate_against(problem)


def test_solve_is_deterministic():
    problem = _problem()
    adapter = MT3DAdapter()
    first = adapter.solve(problem)
    second = adapter.solve(problem)
    np.testing.assert_array_equal(first.impedance_v_a, second.impedance_v_a)


def test_solve_honors_requested_component_subset_and_order():
    problem = _problem(components=("zyx", "zxy"))
    result = MT3DAdapter().solve(problem)
    assert result.components == ("zyx", "zxy")
    assert result.shape == (1, 1, 2)


def test_diagnostics_report_zero_iterations_for_a_direct_solver():
    problem = _problem()
    result = MT3DAdapter().solve(problem)
    assert np.all(result.diagnostics.iterations == 0)
    assert result.diagnostics.runtime_s >= 0.0
    assert np.all(np.isfinite(result.diagnostics.relative_residual))


def test_diagonal_components_vanish_and_off_diagonal_are_antisymmetric():
    # A half-space model: Zxx=Zyy=0 and Zyx=-Zxy, regardless of the
    # boundary-condition approximation's absolute accuracy.
    mesh, receivers = _calibrated_mesh()
    problem = MaxwellProblem(
        mesh,
        np.full(mesh.shape, 1.0 / 100.0),
        [1.0],
        receivers,
        ("zxx", "zxy", "zyx", "zyy"),
    )
    result = MT3DAdapter().solve(problem)
    zxx, zxy, zyx, zyy = result.impedance_v_a[0, 0, :]
    assert abs(zxx) / abs(zxy) < 1e-8
    assert abs(zyy) / abs(zxy) < 1e-8
    np.testing.assert_allclose(zyx, -zxy, rtol=0.05)


# --------------------------------------------------------------------------- #
# M6 gate: what this research prototype actually earns and does not
# --------------------------------------------------------------------------- #


def test_half_space_benchmark_passes_within_default_thresholds():
    mesh, receivers = _calibrated_mesh()
    freqs = [0.3, 1.0, 3.0]
    benchmark = half_space_benchmark(
        mesh, receivers, freqs, resistivity_ohm_m=100.0
    )
    outcome = benchmark.run(MT3DAdapter())
    assert outcome.passed, outcome.failures
    assert outcome.metrics.normalized_rms < 0.05


def test_layered_earth_benchmark_does_not_pass_yet():
    """Documented, tracked limitation (see the module docstring's
    "Measured accuracy" section): the single-layer boundary
    approximation is not accurate enough for a genuinely layered
    earth, and refining the mesh does not fix it. If this test starts
    failing because the benchmark now passes, update this test and the
    module docstring/`_VERIFIED_BENCHMARKS` together -- do not just
    delete the assertion.
    """
    nx, ny, nz = 6, 6, 14
    hx, hy, hz = 200.0, 200.0, 100.0
    mesh = MaxwellMesh(
        np.arange(nx + 1) * hx,
        np.arange(nz + 1) * hz,
        np.arange(ny + 1) * hy,
    )
    receivers = ReceiverSet([[nx * hx / 2, ny * hy / 2, 0.0]], ["S00"])
    z_edges = np.arange(nz + 1) * hz
    interface1, interface2 = z_edges[3], z_edges[7]
    benchmark = layered_earth_benchmark(
        mesh,
        receivers,
        [0.5, 1.0, 2.0],
        [100.0, 30.0, 500.0],
        [interface1, interface2 - interface1],
    )
    outcome = benchmark.run(MT3DAdapter())
    assert not outcome.passed
    assert outcome.metrics.normalized_rms > 0.1
