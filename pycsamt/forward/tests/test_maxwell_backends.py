"""Tests for Maxwell capability declarations and lazy backend registry."""

import numpy as np
import pytest

from pycsamt.forward.maxwell import (
    BackendCapabilities,
    BackendRegistration,
    BackendRegistry,
    ForwardResult,
    MaxwellBackend,
    MaxwellMesh,
    MaxwellProblem,
    ReceiverSet,
    SolverDiagnostics,
    TriMesh,
    TriProblem,
)


def _problem(*, dimension=2, nonuniform=False, inactive=None, components=None):
    x = [0, 1, 3] if nonuniform else [0, 1, 2]
    mesh = (
        MaxwellMesh(x, [0, 1, 2])
        if dimension == 2
        else MaxwellMesh(x, [0, 1, 2], [0, 1, 2])
    )
    receivers = ReceiverSet([[0.5, 0]] if dimension == 2 else [[0.5, 0.5, 0]], ["S"])
    selected = components or (
        ("zxy", "zyx") if dimension == 2 else ("zxx", "zxy", "zyx", "zyy")
    )
    return MaxwellProblem(
        mesh, np.ones(mesh.shape), [10, 1], receivers, selected, inactive
    )


class _DemoBackend:
    def __init__(self, capabilities, scale=1):
        self.capabilities = capabilities
        self.scale = scale

    def solve(self, problem):
        self.capabilities.assess(problem).require()
        diagnostics = SolverDiagnostics(
            np.ones((len(problem.frequencies_hz), 1), bool),
            np.zeros((len(problem.frequencies_hz), 1), int),
            np.zeros((len(problem.frequencies_hz), 1)),
            0,
        )
        values = np.full(
            (
                problem.receivers.count,
                len(problem.frequencies_hz),
                len(problem.components),
            ),
            self.scale + 1j,
        )
        return ForwardResult(
            problem.problem_hash,
            problem.frequencies_hz,
            problem.receivers.names,
            problem.components,
            values,
            None,
            self.capabilities.name,
            self.capabilities.version,
            diagnostics,
        )


def _capabilities(**changes):
    values = dict(
        name="demo",
        version="1",
        dimensions=(2,),
        components=("zxy", "zyx"),
        verified_benchmarks=("half-space",),
    )
    values.update(changes)
    return BackendCapabilities(**values)


def test_capability_roundtrip_and_normalization():
    capabilities = _capabilities(name="Demo_Backend")
    assert capabilities.name == "demo-backend"
    assert BackendCapabilities.from_dict(capabilities.to_dict()) == capabilities
    assert capabilities.supports_component("ZXY")


def test_assess_collects_dimension_component_mesh_and_limit_errors():
    capabilities = _capabilities(
        dimensions=(3,),
        components=("zxx",),
        supports_nonuniform_mesh=False,
        maximum_cells=2,
        maximum_frequencies=1,
    )
    report = capabilities.assess(_problem(nonuniform=True))
    assert not report.compatible
    assert len(report.errors) == 5
    with pytest.raises(ValueError, match="2-D problems are unsupported"):
        report.require()


def test_inactive_cells_and_topography_capabilities_are_distinct():
    active = np.array([[False, True], [True, True]])
    problem = _problem(inactive=active)
    assert "inactive cells" in _capabilities().assess(problem).errors[0]
    report = _capabilities(supports_inactive_cells=True).assess(problem)
    assert "topography" in report.errors[0]
    assert (
        _capabilities(supports_inactive_cells=True, supports_topography=True)
        .assess(problem)
        .compatible
    )


def test_unverified_backend_gets_warning_not_error():
    report = _capabilities(verified_benchmarks=()).assess(_problem())
    assert report.compatible
    assert report.warnings == ("backend declares no verified benchmarks",)


def test_assess_accepts_tri_problem_without_nonuniform_mesh_penalty():
    mesh = TriMesh(
        [[0, 0], [1000, 0], [1000, 500], [0, 500]], [[0, 1, 2], [0, 2, 3]]
    )
    problem = TriProblem(
        mesh, np.ones(mesh.shape), [10, 1], ReceiverSet([[500, 0]], ["S"])
    )
    # supports_nonuniform_mesh=False must not penalize a TriMesh -- it has
    # no cell_widths_m / uniform-vs-nonuniform concept to check.
    report = _capabilities(supports_nonuniform_mesh=False).assess(problem)
    assert report.compatible


def test_assess_rejects_unrecognized_problem_type():
    with pytest.raises(TypeError, match="MaxwellProblem or TriProblem"):
        _capabilities().assess(object())


def test_registration_is_lazy_and_options_reach_factory():
    capabilities = _capabilities()
    calls = []

    def factory(scale=1):
        calls.append(scale)
        return _DemoBackend(capabilities, scale)

    registration = BackendRegistration(capabilities, factory)
    assert calls == []
    backend = registration.create(scale=3)
    assert isinstance(backend, MaxwellBackend)
    assert calls == [3]
    result = backend.solve(_problem())
    assert result.success and np.all(result.impedance_v_a.real == 3)


def test_registration_reports_unavailable_and_factory_contract_errors():
    capabilities = _capabilities()
    unavailable = BackendRegistration(
        capabilities, lambda: None, lambda: (False, "dependency missing")
    )
    assert unavailable.availability() == (False, "dependency missing")
    with pytest.raises(RuntimeError, match="dependency missing"):
        unavailable.create()
    with pytest.raises(RuntimeError, match="violates"):
        BackendRegistration(capabilities, lambda: object()).create()
    with pytest.raises(RuntimeError, match="differ"):
        BackendRegistration(
            capabilities, lambda: _DemoBackend(_capabilities(version="2"))
        ).create()


def test_registry_lifecycle_replace_describe_and_availability():
    registry = BackendRegistry()
    capabilities = _capabilities()
    registration = BackendRegistration(capabilities, lambda: _DemoBackend(capabilities))
    registry.register(registration)
    with pytest.raises(ValueError, match="already registered"):
        registry.register(registration)
    registry.register(registration, replace=True)
    assert registry.names() == ("demo",)
    assert registry.describe()["demo"]["available"]
    assert registry.create("DEMO").capabilities == capabilities
    assert registry.unregister("demo") == registration
    with pytest.raises(KeyError, match="unknown"):
        registry.get("demo")


def test_registry_filters_unavailable_and_handles_probe_failure():
    registry = BackendRegistry()
    cap_a = _capabilities(name="available")
    cap_b = _capabilities(name="broken")
    registry.register(BackendRegistration(cap_a, lambda: _DemoBackend(cap_a)))

    def broken_probe():
        raise ImportError("no solver")

    registry.register(
        BackendRegistration(cap_b, lambda: _DemoBackend(cap_b), broken_probe)
    )
    assert registry.names(available_only=True) == ("available",)
    assert "probe failed" in registry.describe()["broken"]["reason"]
