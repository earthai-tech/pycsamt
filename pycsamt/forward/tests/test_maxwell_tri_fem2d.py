"""Tests for the in-house 2-D triangular FEM MT solver (tri_fem2d.py, M16).

Layered like ``test_maxwell_mt3d.py``: pure linear-algebra/FEM identity
checks first (before any EM physics runs), then preflight, then solve
round-trip sanity, then analytic half-space/layered-earth benchmarks
(checked directly against ``half_space_impedance``/``layered_earth_impedance``,
not through ``MaxwellBenchmark`` -- that class is hard-typed to
``MaxwellMesh`` and building a triangular-specific benchmark class was
deliberately out of scope, same call already made for ``mare2dem.py``),
and finally a cross-check against the independently-validated ``mt2d.py``
finite-difference solver on a shared heterogeneous model.

The benchmark mesh parameters below were chosen empirically, the same way
``test_maxwell_mt2d.py`` documents its own lateral-extent gotcha: accuracy
of the ``Hx``/``Ex`` area-weighted gradient recovery at a receiver depends
on how fine the triangles immediately *below* the receiver are relative to
skin depth, not just overall domain size -- see ``tri_fem2d.py``'s own
module docstring ("Station field extraction") for the numbers. A uniform
mesh at ~0.43 skin depths per cell gave 10-33% error; grading the mesh to
~0.1 skin depths near the surface only (not everywhere -- that would be
needlessly expensive) brought it under 1%.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse.linalg import spsolve
from scipy.spatial import Delaunay

from pycsamt.forward.maxwell import (
    IncompatibleProblemError,
    MaxwellMesh,
    MaxwellProblem,
    ReceiverSet,
    create_backend,
    list_backends,
)
from pycsamt.forward.maxwell.benchmarks import (
    half_space_impedance,
    layered_earth_impedance,
)
from pycsamt.forward.maxwell.contracts_tri import TriMesh, TriProblem
from pycsamt.forward.maxwell.mesh import skin_depth_m
from pycsamt.forward.maxwell.mt2d import MT2DAdapter
from pycsamt.forward.maxwell.tri_fem2d import (
    MU0,
    TriFEM2DAdapter,
    _apply_dirichlet,
    _boundary_node_ids,
    assemble_stiffness_mass,
    p1_gradients,
    register_trifem2d_backend,
)

_MU0 = MU0


# ---------------------------------------------------------------------------
# Mesh-building helpers shared across benchmark tests
# ---------------------------------------------------------------------------


def _graded_z(fine_max: float, depth: float, *, breaks=(), n_fine=21, n_coarse=40):
    zs = np.concatenate(
        [
            np.linspace(0.0, fine_max, n_fine),
            np.linspace(fine_max, depth, n_coarse)[1:],
        ]
    )
    zs = np.round(zs, 6)
    if breaks:
        zs = np.union1d(zs, np.round(np.asarray(breaks, dtype=float), 6))
    return np.unique(zs)


def _uniform_mesh(x_half_width, depth, station_x, *, fine_max=1500.0, nx=25, breaks=()):
    # Rounding before union/unique matters: linspace's own floating-point
    # roundoff can place its middle point a few 1e-12 m away from an
    # explicitly-requested station x (e.g. 0.0), which union1d then keeps
    # as a second, near-duplicate x-column -- Qhull's Delaunay silently
    # produces zero triangles touching that sliver column, leaving those
    # mesh nodes with no incident triangle at all (an exactly-singular row
    # in the assembled system, since they are not boundary nodes either).
    # Found via a real singular-matrix failure while writing this helper.
    xs = np.union1d(
        np.round(np.linspace(-x_half_width, x_half_width, nx), 6),
        np.round(np.asarray(station_x, dtype=float), 6),
    )
    zs = _graded_z(fine_max, depth, breaks=breaks)
    xx, zz = np.meshgrid(xs, zs)
    points = np.unique(np.column_stack([xx.ravel(), zz.ravel()]), axis=0)
    tri = Delaunay(points)
    mesh = TriMesh(points, tri.simplices, boundary_segments=tri.convex_hull)
    return mesh


def _receiver_node_ids(mesh: TriMesh, station_x) -> np.ndarray:
    ids = []
    for x in station_x:
        d = np.hypot(mesh.nodes_m[:, 0] - x, mesh.nodes_m[:, 1])
        ids.append(int(np.argmin(d)))
    return np.array(ids)


def _half_space_setup(resistivity_ohm_m, frequencies_hz, station_x=(0.0,)):
    min_freq = min(frequencies_hz)
    max_freq = max(frequencies_hz)
    half_width = 8.0 * float(skin_depth_m(resistivity_ohm_m, min_freq))
    fine_max = 0.3 * float(skin_depth_m(resistivity_ohm_m, max_freq))
    mesh = _uniform_mesh(half_width, half_width, station_x, fine_max=fine_max)
    node_ids = _receiver_node_ids(mesh, station_x)
    resistivity = np.full(mesh.n_triangles, resistivity_ohm_m)
    return mesh, node_ids, resistivity


def _layered_setup(resistivity_ohm_m, thickness_m, frequencies_hz, station_x=(0.0,)):
    # Lateral padding must scale with the *largest* skin depth present, not
    # the smallest resistivity's -- sizing off the conductive middle layer
    # here first gave a domain too narrow relative to the deep, resistive
    # basal half-space at the lowest frequency, and a 14% TM-mode error at
    # that combination is what exposed it (same class of gotcha
    # test_maxwell_mt2d.py already documents for its own FD solver).
    background_max = max(resistivity_ohm_m)
    min_freq = min(frequencies_hz)
    max_freq = max(frequencies_hz)
    half_width = 8.0 * float(skin_depth_m(background_max, min_freq))
    fine_max = 0.3 * float(skin_depth_m(resistivity_ohm_m[0], max_freq))
    breaks = np.cumsum(thickness_m)
    depth = max(half_width, float(breaks[-1]) * 3.0) if len(breaks) else half_width
    mesh = _uniform_mesh(
        half_width, depth, station_x, fine_max=fine_max, breaks=tuple(breaks)
    )
    node_ids = _receiver_node_ids(mesh, station_x)
    z_top = np.concatenate([[0.0], breaks])
    centroids_z = mesh.triangle_centroids_m[:, 1]
    layer_idx = np.clip(
        np.searchsorted(z_top, centroids_z, side="right") - 1,
        0,
        len(resistivity_ohm_m) - 1,
    )
    resistivity = np.asarray(resistivity_ohm_m)[layer_idx]
    return mesh, node_ids, resistivity


def _problem(mesh: TriMesh, resistivity, frequencies, station_x, *, components=("zxy", "zyx")):
    receivers = ReceiverSet([[x, 0.0] for x in station_x], [f"S{i:02d}" for i in range(len(station_x))])
    return TriProblem(mesh, 1.0 / resistivity, frequencies, receivers, components)


# ---------------------------------------------------------------------------
# Layer 1: pure linear-algebra / FEM identity checks, no physics
# ---------------------------------------------------------------------------


def test_p1_gradients_reproduce_a_known_linear_field_exactly():
    nodes = np.array([[0.0, 0.0], [173.0, -12.0], [40.0, 210.0]])
    triangle = np.array([0, 1, 2])
    grad, area = p1_gradients(nodes, triangle)
    assert area > 0.0

    def f(p):
        return 2.0 * p[0] - 3.0 * p[1] + 7.0

    values = np.array([f(p) for p in nodes])
    recovered_gradient = values @ grad
    np.testing.assert_allclose(recovered_gradient, [2.0, -3.0], atol=1e-10)


def test_p1_gradients_are_orientation_invariant():
    nodes = np.array([[0.0, 0.0], [173.0, -12.0], [40.0, 210.0]])
    grad_ccw, area_ccw = p1_gradients(nodes, np.array([0, 1, 2]))
    grad_cw, area_cw = p1_gradients(nodes, np.array([0, 2, 1]))
    assert area_ccw == pytest.approx(area_cw)
    stiffness_ccw = grad_ccw @ grad_ccw.T
    # reorder the CW triangle's rows/cols back to node order 0,1,2
    perm = [0, 2, 1]
    stiffness_cw = (grad_cw @ grad_cw.T)[np.ix_(perm, perm)]
    np.testing.assert_allclose(stiffness_ccw, stiffness_cw, atol=1e-10)


def test_local_mass_row_sum_equals_area_over_three():
    mesh = TriMesh([[0, 0], [10, 0], [0, 8]], [[0, 1, 2]])
    _, M = assemble_stiffness_mass(
        mesh, np.array([1.0]), weight_stiffness=False
    )
    area = mesh.triangle_areas_m2[0]
    np.testing.assert_allclose(np.asarray(M.sum(axis=1)).ravel(), area / 3.0)


def test_stiffness_matrix_row_sums_are_zero():
    """Constant fields have zero gradient, so K @ ones == 0 exactly."""
    mesh = _uniform_mesh(2000.0, 2000.0, (0.0,), fine_max=500.0, nx=9)
    K, _ = assemble_stiffness_mass(
        mesh, np.ones(mesh.n_triangles), weight_stiffness=False
    )
    np.testing.assert_allclose(
        K @ np.ones(mesh.n_nodes), np.zeros(mesh.n_nodes), atol=1e-8
    )


def test_p1_patch_test_reproduces_linear_dirichlet_data_exactly():
    """A pure Laplace solve with linear Dirichlet data must recover that
    exact linear function at every interior node -- the classic FEM
    correctness check, and the one that caught this module's real sign
    bug during development (see tri_fem2d.py's module docstring).
    """
    mesh = _uniform_mesh(2000.0, 3000.0, (0.0,), fine_max=500.0, nx=11)
    K, _ = assemble_stiffness_mass(
        mesh, np.ones(mesh.n_triangles), weight_stiffness=False
    )

    def f(p):
        return 0.4 * p[0] - 0.7 * p[1] + 3.0

    exact = np.array([f(p) for p in mesh.nodes_m])
    boundary_ids = _boundary_node_ids(mesh)
    b = np.zeros(mesh.n_nodes)
    A, b = _apply_dirichlet(
        K.astype(complex), b.astype(complex), boundary_ids, exact[boundary_ids]
    )
    solved = spsolve(A, b)
    np.testing.assert_allclose(solved.real, exact, atol=1e-6)
    np.testing.assert_allclose(solved.imag, 0.0, atol=1e-9)


def test_vertical_column_profile_reproduces_layered_structure():
    from pycsamt.forward.maxwell.tri_fem2d import _vertical_column_profile

    mesh, _, resistivity = _layered_setup(
        [300.0, 30.0, 1000.0], [500.0, 1500.0], [1.0]
    )
    rho, thick = _vertical_column_profile(mesh, resistivity, 0.0)
    # The recovered layered column should start at the true near-surface
    # resistivity and end at the true basal resistivity.
    assert rho[0] == pytest.approx(300.0)
    assert rho[-1] == pytest.approx(1000.0)
    assert thick.sum() > 1900.0  # covers at least the two finite layers


# ---------------------------------------------------------------------------
# Layer 2: assess() preflight
# ---------------------------------------------------------------------------


def _small_mesh(*, with_boundary=True):
    return TriMesh(
        [[0, 0], [1000, 0], [1000, 500], [0, 500]],
        [[0, 1, 2], [0, 2, 3]],
        boundary_segments=(
            [[0, 1], [1, 2], [2, 3], [3, 0]] if with_boundary else None
        ),
    )


def test_assess_accepts_a_valid_surface_receiver_problem():
    mesh = _small_mesh()
    problem = TriProblem(
        mesh, np.full(2, 0.01), [1.0], ReceiverSet([[0.0, 0.0]], ["S00"])
    )
    report = TriFEM2DAdapter().assess(problem)
    assert report.compatible


def test_assess_rejects_missing_boundary_segments():
    mesh = _small_mesh(with_boundary=False)
    problem = TriProblem(
        mesh, np.full(2, 0.01), [1.0], ReceiverSet([[0.0, 0.0]], ["S00"])
    )
    report = TriFEM2DAdapter().assess(problem)
    assert not report.compatible
    assert "boundary_segments" in report.errors[0]


def test_assess_rejects_receiver_not_on_a_mesh_node():
    mesh = _small_mesh()
    problem = TriProblem(
        mesh, np.full(2, 0.01), [1.0], ReceiverSet([[500.0, 0.0]], ["S00"])
    )
    report = TriFEM2DAdapter().assess(problem)
    assert not report.compatible
    assert "coincide with a mesh node" in report.errors[-1]


def test_assess_rejects_buried_receivers():
    mesh = _small_mesh()
    problem = TriProblem(
        mesh, np.full(2, 0.01), [1.0], ReceiverSet([[1000.0, 500.0]], ["S00"])
    )
    report = TriFEM2DAdapter().assess(problem)
    assert not report.compatible
    assert any("surface" in e for e in report.errors)


def test_assess_rejects_non_vacuum_permeability():
    mesh = _small_mesh()
    problem = TriProblem(
        mesh,
        np.full(2, 0.01),
        [1.0],
        ReceiverSet([[0.0, 0.0]], ["S00"]),
        magnetic_permeability_h_m=1.0,
    )
    report = TriFEM2DAdapter().assess(problem)
    assert not report.compatible
    assert any("permeability" in e for e in report.errors)


def test_assess_rejects_a_rectilinear_maxwell_problem():
    mesh = MaxwellMesh([0, 500, 1000], [0, 250, 500])
    problem = MaxwellProblem(
        mesh, np.full(mesh.shape, 0.01), [1.0], ReceiverSet([[500.0, 0.0]], ["S00"])
    )
    report = TriFEM2DAdapter().assess(problem)
    assert not report.compatible


def test_solve_raises_incompatible_before_writing_any_file():
    mesh = _small_mesh(with_boundary=False)
    problem = TriProblem(
        mesh, np.full(2, 0.01), [1.0], ReceiverSet([[0.0, 0.0]], ["S00"])
    )
    with pytest.raises(IncompatibleProblemError, match="boundary_segments"):
        TriFEM2DAdapter().solve(problem)


# ---------------------------------------------------------------------------
# Layer 3: solve round-trip sanity
# ---------------------------------------------------------------------------


def test_solve_returns_canonical_shape_and_is_deterministic():
    mesh, node_ids, resistivity = _half_space_setup(100.0, [1.0, 0.5], (0.0,))
    problem = _problem(mesh, resistivity, [1.0, 0.5], (0.0,))
    adapter = TriFEM2DAdapter()
    result_a = adapter.solve(problem)
    result_b = adapter.solve(problem)
    assert result_a.shape == (1, 2, 2)
    assert result_a.success
    assert result_a.backend_name == "trifem2d"
    np.testing.assert_array_equal(
        result_a.impedance_v_a, result_b.impedance_v_a
    )


def test_solve_honors_component_subset_and_order():
    mesh, node_ids, resistivity = _half_space_setup(100.0, [1.0], (0.0,))
    problem = _problem(mesh, resistivity, [1.0], (0.0,), components=("zyx",))
    result = TriFEM2DAdapter().solve(problem)
    assert result.components == ("zyx",)
    assert result.shape == (1, 1, 1)


def test_solve_diagnostics_are_finite_and_direct():
    mesh, node_ids, resistivity = _half_space_setup(100.0, [1.0], (0.0,))
    problem = _problem(mesh, resistivity, [1.0], (0.0,))
    result = TriFEM2DAdapter().solve(problem)
    assert np.all(result.diagnostics.converged)
    assert np.all(result.diagnostics.iterations == 0)
    assert np.all(np.isfinite(result.diagnostics.relative_residual))
    assert result.diagnostics.runtime_s >= 0.0


def test_solve_half_space_gives_zyx_close_to_minus_zxy():
    mesh, node_ids, resistivity = _half_space_setup(100.0, [1.0], (0.0,))
    problem = _problem(mesh, resistivity, [1.0], (0.0,))
    result = TriFEM2DAdapter().solve(problem)
    zxy = result.impedance_v_a[0, 0, 0]
    zyx = result.impedance_v_a[0, 0, 1]
    np.testing.assert_allclose(zyx, -zxy, rtol=0.05)


# ---------------------------------------------------------------------------
# Layer 4: analytic benchmarks
# ---------------------------------------------------------------------------


def test_half_space_benchmark_matches_analytic_impedance():
    resistivity_ohm_m = 100.0
    frequencies = [1.0, 0.1]
    station_x = (-500.0, 0.0, 500.0)
    mesh, node_ids, resistivity = _half_space_setup(
        resistivity_ohm_m, frequencies, station_x
    )
    problem = _problem(mesh, resistivity, frequencies, station_x)
    result = TriFEM2DAdapter().solve(problem)

    errors = []
    for fi, freq in enumerate(frequencies):
        expected_zxy = half_space_impedance(resistivity_ohm_m, freq)
        for ri in range(len(station_x)):
            zxy = result.impedance_v_a[ri, fi, 0]
            zyx = result.impedance_v_a[ri, fi, 1]
            errors.append(abs(zxy - expected_zxy) / abs(expected_zxy))
            errors.append(abs(zyx - (-expected_zxy)) / abs(expected_zxy))
    assert max(errors) < 0.05, f"max relative error {max(errors):.4f}"


def test_half_space_benchmark_on_a_real_graded_triangle_mesh():
    """Same half-space benchmark as above, but on a mesh built by
    :func:`~pycsamt.forward.maxwell.tri_mesh_gen.build_graded_tri_mesh`
    (real Shewchuk-Triangle quality/grading) instead of this file's own
    hand-built structured mesh -- confirms the mesh-generation redesign
    (replacing ``dataset2d_tri.py``'s old uniform-lattice
    ``scipy.spatial.Delaunay`` approach) still gives correct FEM physics,
    not just a better-looking picture. ``surface_cell_m`` must be chosen
    relative to skin depth for solver accuracy, exactly as
    ``tri_fem2d.py``'s own "Station field extraction" docstring section
    describes for any triangular mesh, graded or structured.
    """
    from pycsamt.forward.maxwell.tri_mesh_gen import build_graded_tri_mesh

    resistivity_ohm_m = 100.0
    frequencies = [1.0, 0.1]
    station_x = (-500.0, 0.0, 500.0)
    half_width = 8.0 * float(skin_depth_m(resistivity_ohm_m, min(frequencies)))
    surface_cell_m = 0.04 * float(skin_depth_m(resistivity_ohm_m, max(frequencies)))

    mesh = build_graded_tri_mesh(
        (-half_width, half_width),
        (0.0, half_width),
        station_x,
        surface_cell_m=surface_cell_m,
        growth_rate=1.3,
    )
    resistivity = np.full(mesh.n_triangles, resistivity_ohm_m)
    problem = _problem(mesh, resistivity, frequencies, station_x)
    result = TriFEM2DAdapter().solve(problem)

    errors = []
    for fi, freq in enumerate(frequencies):
        expected_zxy = half_space_impedance(resistivity_ohm_m, freq)
        for ri in range(len(station_x)):
            zxy = result.impedance_v_a[ri, fi, 0]
            zyx = result.impedance_v_a[ri, fi, 1]
            errors.append(abs(zxy - expected_zxy) / abs(expected_zxy))
            errors.append(abs(zyx - (-expected_zxy)) / abs(expected_zxy))
    assert max(errors) < 0.05, f"max relative error {max(errors):.4f}"


def test_half_space_benchmark_is_invariant_to_a_constant_vertical_datum_shift():
    """Maxwell's equations don't care about an absolute vertical datum --
    translating the whole mesh + receivers by a constant offset must
    reproduce the identical impedance. This isolates and proves the
    ``tri_fem2d.py`` surface-detection fix (``_local_surface_z``,
    replacing the old hardcoded ``abs(z)<=tol`` assumption) without
    needing an independent topographic-MT reference solution -- see
    ``tri_fem2d.py``'s own "Scope, honestly stated" section for why
    genuine spatially-varying topographic distortion itself is not
    benchmarked here.
    """
    from pycsamt.forward.maxwell.tri_mesh_gen import build_graded_tri_mesh

    resistivity_ohm_m = 100.0
    frequencies = [1.0]
    station_x = (-500.0, 0.0, 500.0)
    half_width = 8.0 * float(skin_depth_m(resistivity_ohm_m, min(frequencies)))
    surface_cell_m = 0.04 * float(skin_depth_m(resistivity_ohm_m, max(frequencies)))

    def _max_relative_error(z_shift: float) -> float:
        mesh = build_graded_tri_mesh(
            (-half_width, half_width),
            (0.0, half_width),
            station_x,
            surface_cell_m=surface_cell_m,
            growth_rate=1.3,
        )
        nodes = mesh.nodes_m.copy()
        nodes[:, 1] += z_shift
        shifted_mesh = TriMesh(
            nodes,
            mesh.triangles,
            region_ids=mesh.region_ids,
            boundary_segments=mesh.boundary_segments,
        )
        resistivity = np.full(shifted_mesh.n_triangles, resistivity_ohm_m)
        receivers = ReceiverSet(
            [[x, z_shift] for x in station_x], ["S0", "S1", "S2"]
        )
        problem = TriProblem(shifted_mesh, 1.0 / resistivity, frequencies, receivers)
        result = TriFEM2DAdapter().solve(problem)
        expected = half_space_impedance(resistivity_ohm_m, frequencies[0])
        return max(
            abs(result.impedance_v_a[i, 0, 0] - expected) / abs(expected)
            for i in range(len(station_x))
        )

    baseline_error = _max_relative_error(0.0)
    shifted_error = _max_relative_error(200.0)
    assert baseline_error < 0.05
    assert shifted_error == pytest.approx(baseline_error, abs=1e-8)


def test_assess_accepts_receiver_on_a_non_flat_local_surface():
    """A receiver sitting at the mesh's own local top elevation (not
    ``z=0``) must be accepted, and one that is genuinely off that
    surface must still be rejected -- proving the fix generalizes rather
    than merely disabling the check.
    """
    # A simple two-triangle mesh whose "surface" (top edge) sits at
    # z=-50 (elevated terrain, not the origin) rather than z=0.
    mesh = TriMesh(
        [[0, -50], [1000, -50], [1000, 500], [0, 500]],
        [[0, 1, 2], [0, 2, 3]],
        boundary_segments=[[0, 1], [1, 2], [2, 3], [3, 0]],
    )
    on_surface = TriProblem(
        mesh, np.full(2, 0.01), [1.0], ReceiverSet([[0.0, -50.0]], ["S00"])
    )
    report = TriFEM2DAdapter().assess(on_surface)
    assert report.compatible

    # (1000, 500) is a real mesh node (the bottom-right corner) but not on
    # the local surface (which sits at z=-50 for every x here) -- isolates
    # the surface check from the separate node-matching check.
    off_surface = TriProblem(
        mesh, np.full(2, 0.01), [1.0], ReceiverSet([[1000.0, 500.0]], ["S00"])
    )
    report = TriFEM2DAdapter().assess(off_surface)
    assert not report.compatible
    assert any("local surface" in e for e in report.errors)


def test_layered_earth_benchmark_matches_analytic_impedance():
    resistivity_ohm_m = [300.0, 30.0, 1000.0]
    thickness_m = [500.0, 1500.0]
    frequencies = [1.0, 0.1]
    station_x = (0.0,)
    mesh, node_ids, resistivity = _layered_setup(
        resistivity_ohm_m, thickness_m, frequencies, station_x
    )
    problem = _problem(mesh, resistivity, frequencies, station_x)
    result = TriFEM2DAdapter().solve(problem)

    errors = []
    for fi, freq in enumerate(frequencies):
        expected_zxy = layered_earth_impedance(
            resistivity_ohm_m, thickness_m, [freq]
        )[0]
        zxy = result.impedance_v_a[0, fi, 0]
        zyx = result.impedance_v_a[0, fi, 1]
        errors.append(abs(zxy - expected_zxy) / abs(expected_zxy))
        errors.append(abs(zyx - (-expected_zxy)) / abs(expected_zxy))
    assert max(errors) < 0.07, f"max relative error {max(errors):.4f}"


# ---------------------------------------------------------------------------
# Layer 5: cross-check against the independently-validated mt2d.py
# ---------------------------------------------------------------------------


def test_cross_check_against_mt2d_on_a_lateral_anomaly():
    """Both solvers see the same background + conductive block; their
    surface responses at shared stations should agree to a loose (cross-
    check, not precision-benchmark) tolerance -- the validation matrix's
    own "backend cross-check" category, and the only ground truth this
    module has for laterally-heterogeneous structure.
    """
    background = 100.0
    anomaly = 10.0
    freq = 1.0
    station_x = (-1000.0, 0.0, 1000.0)
    half_width = 8.0 * float(skin_depth_m(background, freq))

    # -- mt2d.py (rectilinear FD) --
    # Grid2D (which mt2d.py's adapter builds from mesh.cell_widths_m) always
    # starts its internal domain at x=0, discarding a MaxwellMesh's actual
    # edge coordinates -- confirmed empirically while writing this test (a
    # centered [-half_width, +half_width] mesh solved receivers against an
    # internally-shifted [0, 2*half_width] grid, "off the grid" for any
    # negative station x). Shift this solver's own coordinate frame to
    # [0, 2*half_width] to match that convention; the triangular mesh below
    # keeps the natural centered frame, since TriMesh has no such quirk.
    dx = np.full(48, 2.0 * half_width / 48)
    x_edges = np.concatenate([[0.0], np.cumsum(dx)])
    dz = np.full(40, half_width / 40)
    z_edges = np.concatenate([[0.0], np.cumsum(dz)])
    rect_mesh = MaxwellMesh(x_edges, z_edges)
    conductivity = np.full(rect_mesh.shape, 1.0 / background)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:]) - half_width
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    anomaly_mask = (
        (np.abs(x_centers) < 300.0)[None, :]
        & (z_centers < 800.0)[:, None]
    )
    conductivity[anomaly_mask] = 1.0 / anomaly
    rect_problem = MaxwellProblem(
        rect_mesh,
        conductivity,
        [freq],
        ReceiverSet(
            [[x + half_width, 0.0] for x in station_x], ["S0", "S1", "S2"]
        ),
    )
    rect_result = MT2DAdapter(verbose=False).solve(rect_problem)

    # -- trifem2d (triangular FEM), same background + anomaly footprint --
    fine_max = 0.3 * float(skin_depth_m(background, freq))
    mesh = _uniform_mesh(half_width, half_width, station_x, fine_max=fine_max)
    centroids = mesh.triangle_centroids_m
    resistivity = np.full(mesh.n_triangles, background)
    tri_anomaly_mask = (np.abs(centroids[:, 0]) < 300.0) & (
        centroids[:, 1] < 800.0
    )
    resistivity[tri_anomaly_mask] = anomaly
    problem = _problem(mesh, resistivity, [freq], station_x)
    tri_result = TriFEM2DAdapter().solve(problem)

    for ri in range(len(station_x)):
        for ci in range(2):
            rect_z = rect_result.impedance_v_a[ri, 0, ci]
            tri_z = tri_result.impedance_v_a[ri, 0, ci]
            # loose cross-check tolerance: independent discretizations
            # (FD tensor grid vs. unstructured triangles) of the same
            # sharp-edged block will not agree to benchmark precision.
            assert abs(tri_z - rect_z) / abs(rect_z) < 0.35, (
                ri, ci, rect_z, tri_z
            )


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------


def test_register_trifem2d_backend_is_immediately_available():
    register_trifem2d_backend(replace=True)
    assert "trifem2d" in list_backends()
    assert list_backends()["trifem2d"]["available"] is True
    backend = create_backend("trifem2d")
    assert backend.capabilities.name == "trifem2d"
