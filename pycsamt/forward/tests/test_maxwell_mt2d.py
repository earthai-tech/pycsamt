"""Tests for the validated 2-D MT finite-difference adapter (mt2d.py).

The benchmark mesh parameters below were chosen empirically: this
solver applies Dirichlet boundary conditions at the lateral edges from
the exact 1-D response of the edge columns, so accuracy on a laterally
uniform model depends on the *lateral* domain extent relative to skin
depth, not only on depth/z resolution. A domain that is deep enough
but too narrow left 25-45% amplitude errors at the lowest frequency
during calibration; widening x to comfortably exceed the deepest skin
depth (not just z) is what made both analytic benchmarks pass.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    BackendRegistry,
    CompatibilityReport,
    IncompatibleProblemError,
    MaxwellMesh,
    MaxwellProblem,
    ReceiverSet,
    create_backend,
    half_space_benchmark,
    layered_earth_benchmark,
    list_backends,
    run_benchmarks,
)
from pycsamt.forward.maxwell.backends import BackendRegistration
from pycsamt.forward.maxwell.mt2d import MT2DAdapter, register_mt2d_backend

_MU0 = 4.0e-7 * np.pi


def _calibrated_mesh():
    """Return the (mesh, receivers) pair calibrated for both benchmarks."""
    dz0, growth, n = 25.0, 1.22, 34
    dz = dz0 * growth ** np.arange(n)
    z_edges = np.concatenate([[0.0], np.cumsum(dz)])
    x_edges = np.linspace(0, 240_000, 25)
    mesh = MaxwellMesh(x_edges, z_edges)
    receivers = ReceiverSet([[120_000.0, 0.0]], ["S00"])
    return mesh, receivers, z_edges


def _problem(**changes):
    mesh = MaxwellMesh([0, 5_000, 10_000], [0, 2_500, 5_000])
    values = dict(
        mesh=mesh,
        conductivity_s_m=np.full(mesh.shape, 0.01),
        frequencies_hz=[10.0, 1.0],
        receivers=ReceiverSet([[5_000.0, 0.0]], ["S00"]),
        components=("zxy", "zyx"),
    )
    values.update(changes)
    return MaxwellProblem(**values)


# --------------------------------------------------------------------------- #
# Capability declaration and assess() preflight
# --------------------------------------------------------------------------- #


def test_capabilities_declare_2d_te_tm_only():
    adapter = MT2DAdapter()
    cap = adapter.capabilities
    assert cap.name == "mt2d"
    assert cap.dimensions == (2,)
    assert set(cap.components) == {"zxy", "zyx"}
    assert cap.time_conventions == ("exp(+iwt)",)
    assert cap.supports_inactive_cells is False
    assert cap.supports_topography is False


def test_assess_rejects_buried_receivers():
    problem = _problem(
        receivers=ReceiverSet([[5_000.0, 100.0]], ["S00"])
    )
    report = MT2DAdapter().assess(problem)
    assert not report.compatible
    assert any("surface" in message for message in report.errors)


def test_assess_rejects_receivers_outside_mesh_x_bounds():
    mesh = MaxwellMesh([0, 5_000, 10_000], [0, 2_500, 5_000])
    problem = _problem(
        mesh=mesh,
        receivers=ReceiverSet([[50_000.0, 0.0]], ["S00"]),
    )
    report = MT2DAdapter().assess(problem)
    assert not report.compatible
    assert any("mesh" in message for message in report.errors)


def test_assess_rejects_non_vacuum_permeability():
    problem = _problem(magnetic_permeability_h_m=_MU0 * 2.0)
    report = MT2DAdapter().assess(problem)
    assert not report.compatible
    assert any("permeability" in message for message in report.errors)


def test_assess_combines_base_and_solver_specific_errors():
    problem = _problem(
        receivers=ReceiverSet([[5_000.0, 50.0]], ["S00"]),
        time_dependence="exp(-iwt)",
    )
    report = MT2DAdapter().assess(problem)
    assert isinstance(report, CompatibilityReport)
    assert not report.compatible
    assert len(report.errors) >= 2


def test_assess_returns_the_base_report_unchanged_when_compatible():
    problem = _problem()
    adapter = MT2DAdapter()
    report = adapter.assess(problem)
    assert report.compatible
    assert report == adapter.capabilities.assess(problem)


def test_incompatible_3d_problem_is_rejected_before_solving():
    mesh = MaxwellMesh([0, 1, 2], [0, 1, 2], [0, 1, 2])
    problem = MaxwellProblem(
        mesh,
        np.full(mesh.shape, 0.01),
        [1.0],
        ReceiverSet([[0.5, 0.5, 0.0]], ["S"]),
        ("zxy",),
    )
    with pytest.raises(IncompatibleProblemError, match="3-D"):
        MT2DAdapter().solve(problem)


def test_inactive_cells_are_rejected():
    mesh = MaxwellMesh([0, 5_000, 10_000], [0, 2_500, 5_000])
    active = np.ones(mesh.shape, dtype=bool)
    active[0, 0] = False
    problem = _problem(mesh=mesh, active_cells=active)
    with pytest.raises(IncompatibleProblemError, match="inactive"):
        MT2DAdapter().solve(problem)


# --------------------------------------------------------------------------- #
# Solve round-trip
# --------------------------------------------------------------------------- #


def test_solve_returns_canonical_shape_and_succeeds():
    problem = _problem()
    result = MT2DAdapter().solve(problem)
    assert result.shape == (1, 2, 2)
    assert result.success
    assert result.backend_name == "mt2d"
    result.validate_against(problem)


def test_solve_is_deterministic():
    problem = _problem()
    adapter = MT2DAdapter()
    first = adapter.solve(problem)
    second = adapter.solve(problem)
    np.testing.assert_array_equal(
        first.impedance_v_a, second.impedance_v_a
    )


def test_solve_honors_requested_component_subset_and_order():
    problem = _problem(components=("zyx",))
    result = MT2DAdapter().solve(problem)
    assert result.components == ("zyx",)
    assert result.shape == (1, 2, 1)


def test_diagnostics_report_zero_iterations_for_a_direct_solver():
    problem = _problem()
    result = MT2DAdapter().solve(problem)
    assert np.all(result.diagnostics.iterations == 0)
    assert result.diagnostics.runtime_s >= 0.0
    assert np.all(np.isfinite(result.diagnostics.relative_residual))


def test_te_and_tm_components_are_opposite_signed_for_a_half_space():
    # Zyx = -Zxy for a 1-D earth; regression guard for the sign bug
    # fixed in pycsamt.forward.em2d._surface_impedance_tm. Uses the
    # calibrated mesh (see module docstring): the small toy fixture in
    # _problem() is too coarse/narrow for this equality to hold.
    mesh, receivers, _ = _calibrated_mesh()
    problem = MaxwellProblem(
        mesh,
        np.full(mesh.shape, 1.0 / 100.0),
        [1.0],
        receivers,
        ("zxy", "zyx"),
    )
    result = MT2DAdapter().solve(problem)
    zxy = result.impedance_v_a[0, :, 0]
    zyx = result.impedance_v_a[0, :, 1]
    np.testing.assert_allclose(zxy, -zyx, rtol=0.05)


# --------------------------------------------------------------------------- #
# Registry integration
# --------------------------------------------------------------------------- #


def test_register_mt2d_backend_is_usable_through_the_registry():
    registry = BackendRegistry()
    registry.register(
        BackendRegistration(MT2DAdapter().capabilities, MT2DAdapter)
    )
    backend = registry.create("mt2d")
    result = backend.solve(_problem())
    assert result.success


def test_register_mt2d_backend_process_wide_registry():
    register_mt2d_backend(replace=True)
    assert "mt2d" in list_backends()
    backend = create_backend("mt2d")
    assert isinstance(backend, MT2DAdapter)
    assert backend.solve(_problem()).success


# --------------------------------------------------------------------------- #
# M4 gate: analytic half-space / layered-earth benchmarks
# --------------------------------------------------------------------------- #


def test_half_space_benchmark_passes_within_default_thresholds():
    mesh, receivers, _ = _calibrated_mesh()
    freqs = np.logspace(-1, 1, 6)
    benchmark = half_space_benchmark(
        mesh, receivers, freqs, resistivity_ohm_m=100.0
    )
    outcome = benchmark.run(MT2DAdapter())
    assert outcome.passed, outcome.failures
    assert outcome.metrics.normalized_rms < 0.05


def test_layered_earth_benchmark_passes_within_default_thresholds():
    mesh, receivers, z_edges = _calibrated_mesh()
    freqs = np.logspace(-1, 1, 6)
    interface1, interface2 = z_edges[8], z_edges[16]
    benchmark = layered_earth_benchmark(
        mesh,
        receivers,
        freqs,
        [100.0, 30.0, 500.0],
        [interface1, interface2 - interface1],
    )
    outcome = benchmark.run(MT2DAdapter())
    assert outcome.passed, outcome.failures
    assert outcome.metrics.normalized_rms < 0.05


def test_run_benchmarks_reports_a_fully_passing_suite():
    mesh, receivers, z_edges = _calibrated_mesh()
    freqs = np.logspace(-1, 1, 6)
    interface1, interface2 = z_edges[8], z_edges[16]
    cases = (
        half_space_benchmark(mesh, receivers, freqs, resistivity_ohm_m=100.0),
        layered_earth_benchmark(
            mesh,
            receivers,
            freqs,
            [100.0, 30.0, 500.0],
            [interface1, interface2 - interface1],
        ),
    )
    report = run_benchmarks(MT2DAdapter(), cases)
    assert report.passed
    assert report.pass_fraction == 1.0
