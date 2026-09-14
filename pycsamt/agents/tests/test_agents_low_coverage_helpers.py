"""Unit coverage for pure helpers used by the heavier agent workflows."""

from __future__ import annotations

from types import SimpleNamespace

import matplotlib
import numpy as np

matplotlib.use("Agg")

from pycsamt.agents._base import AgentResult


class _Layered:
    def __init__(self, resistivity, thickness):
        self.resistivity = resistivity
        self.thickness = thickness


class _Grid2D:
    @classmethod
    def halfspace(cls, **kwargs):
        return ("halfspace", kwargs)

    @classmethod
    def with_anomaly(cls, **kwargs):
        return ("anomaly", kwargs)

    @classmethod
    def from_1d_layers(cls, model, **kwargs):
        return ("layers", model, kwargs)


class _Grid3D:
    @classmethod
    def halfspace(cls, **kwargs):
        return ("halfspace", kwargs)

    @classmethod
    def block_anomaly(cls, **kwargs):
        return ("block", kwargs)


def test_forward_grid_builders_cover_all_model_types():
    from pycsamt.agents.forward import (
        _build_grid_2d, _build_grid_3d, _build_layered_model,
    )

    warnings = []
    model = _build_layered_model(None, _Layered, warnings)
    assert isinstance(model, _Layered) and warnings
    assert _build_layered_model(model, _Layered, []) is model

    supplied = object()
    assert _build_grid_2d({"grid": supplied}, _Layered, _Grid2D, [])[0] is supplied
    assert _build_grid_2d({"model": {"type": "halfspace"}}, _Layered, _Grid2D, [])[0][0] == "halfspace"
    assert _build_grid_2d({"model": {"type": "anomaly"}}, _Layered, _Grid2D, [])[0][0] == "anomaly"
    layers, layered = _build_grid_2d(
        {"model": {"resistivities": [10, 100], "thicknesses": [50]}},
        _Layered, _Grid2D, [],
    )
    assert layers[0] == "layers" and isinstance(layered, _Layered)
    warnings = []
    assert _build_grid_2d({"model": {"type": "odd"}}, _Layered, _Grid2D, warnings)[0][0] == "halfspace"
    assert warnings

    assert _build_grid_3d({"grid": supplied}, _Grid3D, []) is supplied
    assert _build_grid_3d({}, _Grid3D, [])[0] == "halfspace"
    assert _build_grid_3d({"model": {"type": "block"}}, _Grid3D, [])[0] == "block"
    warnings = []
    assert _build_grid_3d({"model": {"type": "odd"}}, _Grid3D, warnings)[0] == "halfspace"
    assert warnings


def test_forward_builder_failures_and_misc(monkeypatch):
    from pycsamt.agents.forward import (
        _build_grid_2d, _build_grid_3d, _build_layered_model, _unwrap_fig,
    )

    class Broken:
        def __init__(self, **kwargs):
            raise ValueError("bad")

        @classmethod
        def halfspace(cls, **kwargs):
            raise ValueError("bad grid")

    assert isinstance(_build_layered_model({}, Broken, []), AgentResult)
    assert isinstance(_build_grid_2d({}, _Layered, Broken, [])[0], AgentResult)
    assert isinstance(_build_grid_3d({}, Broken, []), AgentResult)
    assert _unwrap_fig(None) is None
    fig = SimpleNamespace(savefig=lambda: None)
    assert _unwrap_fig(fig) is fig
    axes = SimpleNamespace(get_figure=lambda: fig)
    assert _unwrap_fig(axes) is fig
    assert _unwrap_fig(object()) is None


def test_forward_rms_helper(monkeypatch):
    from pycsamt.agents.forward import _compute_rms_1d
    import pycsamt.emtools._core as core

    z = np.ones((3, 2, 2), dtype=complex)
    ed = SimpleNamespace(rho=np.full((3, 2, 2), 10.0))
    monkeypatch.setattr(core, "ensure_sites", lambda value, verbose=0: [ed])
    monkeypatch.setattr(core, "_iter_items", iter)
    monkeypatch.setattr(core, "_get_z_block", lambda item: (None, z, np.array([1., 10., 100.])))
    response = SimpleNamespace(rho_a=np.full((3, 2, 2), 10.0))
    assert _compute_rms_1d(object(), response, np.array([1., 10., 100.]), 0, 1) < 1e-12


def test_denoising_helpers(monkeypatch):
    import pycsamt.agents.denoising as mod
    import pycsamt.emtools._core as core
    import pycsamt.emtools.remove_noise as noise

    z = np.ones((4, 2, 2), complex)
    monkeypatch.setattr(core, "_iter_items", iter)
    monkeypatch.setattr(core, "_get_z_block", lambda ed: (None, z if ed else None, np.arange(4)))
    assert mod._compute_snr_proxy([1, 0]).size == 4
    monkeypatch.setattr(noise, "rpca_offdiag_denoise", lambda *a, **k: None)
    monkeypatch.setattr(noise, "hampel_filter_freq", lambda *a, **k: "filtered")
    sites = object()
    assert mod._apply_rpca(sites, 2, []) is sites
    assert mod._apply_hampel(sites, 2, []) == "filtered"
    monkeypatch.setattr(noise, "rpca_offdiag_denoise", lambda *a, **k: 1 / 0)
    warnings = []
    assert mod._apply_rpca(sites, 2, warnings) is sites and warnings


def test_edi_per_item_export_variants(tmp_path, monkeypatch):
    from pycsamt.agents.edi_export import _per_item_export
    import pycsamt.emtools._core as core

    class NewWriter:
        def write_new_edi(self, **kwargs):
            return "new.edi"

    class Writer:
        def write(self, **kwargs):
            return None

    class Wrapper:
        def to_edi(self):
            return Writer()

    class Broken:
        def write(self, **kwargs):
            raise OSError("disk")

    items = [NewWriter(), Wrapper(), object(), Broken()]
    monkeypatch.setattr(core, "_iter_items", lambda sites: iter(items))
    monkeypatch.setattr(core, "_name", lambda ed, i: f"S{i}")
    warnings = []
    written, failed = _per_item_export(items, str(tmp_path), "{station}.edi", True, warnings)
    assert len(written) == 2 and len(failed) == 2 and warnings
    existing = tmp_path / "S0.edi"
    existing.touch()
    monkeypatch.setattr(core, "_iter_items", lambda sites: iter([NewWriter()]))
    written, failed = _per_item_export([], str(tmp_path), "{station}.edi", False, warnings)
    assert not written and failed[0][1] == "file exists"


def test_frequency_selection_plot(monkeypatch):
    from pycsamt.agents.freq_decimation import _plot_selection_summary
    import pycsamt.emtools._core as core

    assert _plot_selection_summary({}, {}, [], 3, 2.0) is None
    ed = object()
    monkeypatch.setattr(core, "_iter_items", lambda sites: iter([ed]))
    monkeypatch.setattr(core, "_name", lambda item, i: "S1")
    monkeypatch.setattr(core, "_get_z_block", lambda item: (None, None, np.array([1., 10., 100.])))
    fig = _plot_selection_summary(
        {"S1": np.array([0.1, 1.0])}, {"S1": np.array([False, True, False])},
        [ed], 3, 2.0,
    )
    assert fig is not None


def test_interpretation_object_llm_and_boundaries(monkeypatch):
    from pycsamt.agents.interpretation import InterpretationAgent, resistivity_to_lithology

    assert resistivity_to_lithology(-1) == "unknown lithology"
    assert "brine" in resistivity_to_lithology(1)
    assert "basement" in resistivity_to_lithology(1e8)
    assert InterpretationAgent().execute({"model": {}}).status == "failed"
    model = SimpleNamespace(resistivity=[5, 100, 2000], thickness=[20, 30])
    agent = InterpretationAgent(api_key="key")
    monkeypatch.setattr(agent, "query_llm", lambda *a, **k: "geological answer")
    result = agent.execute({"layered_model": model, "rms": 1.2, "context": "test"})
    assert result.status == "success" and result.llm_interpretation == "geological answer"


def test_inversion_backend_mocked_success_and_fallback(tmp_path, monkeypatch):
    from pycsamt.agents.inversion_backend import InversionBackendAgent, _plot_inversion_section
    import pycsamt.emtools._core as core
    import pycsamt.inversion as inversion

    monkeypatch.setattr(core, "ensure_sites", lambda value, verbose=0: value)
    monkeypatch.setattr(inversion, "available_backends", lambda: ["builtin"])

    class Config:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    inv_result = SimpleNamespace(
        rms=0.2, n_iter=3, model={"log_rho_section": [1., 2.]},
        station_names=["S1"], history=None,
    )
    monkeypatch.setattr(inversion, "InversionConfig", Config)
    monkeypatch.setattr(inversion, "run_inversion", lambda cfg: inv_result)
    agent = InversionBackendAgent(api_key="key", backend="missing")
    monkeypatch.setattr(agent, "query_llm", lambda *a, **k: "good fit")
    monkeypatch.setattr(agent, "_save_figure", lambda *a, **k: str(tmp_path / "fig.png"))
    result = agent.execute({"sites": object(), "output_dir": str(tmp_path)})
    assert result.status == "success" and result.data["backend"] == "builtin"
    assert result.data["log_rho_section"].shape == (2, 1)
    assert result.llm_interpretation == "good fit" and result.warnings
    assert _plot_inversion_section(np.empty((2, 0)), [], np.nan, "x", "1d") is None


def test_inv3d_small_helpers():
    from pycsamt.agents.inv3d_agent import (
        _agent_thicknesses, _extract_station_xy, _is_profile_geometry,
        _nearest_index, _pad_or_trim, _profile_distance_km, _smooth_depth_axis,
    )

    assert _agent_thicknesses(1, np.array([1.]), None).size == 0
    assert np.array_equal(_pad_or_trim(np.ones((2, 3)), 2), np.ones((2, 2)))
    assert _pad_or_trim(np.ones((2, 2)), 4).shape == (2, 4)
    assert _nearest_index(np.array([0, 5, 10]), 7) == 1
    assert np.allclose(_extract_station_xy(SimpleNamespace(), 0), [0, 0])
    assert np.all(np.isfinite(_extract_station_xy(SimpleNamespace(lat=2, lon=3), 0)))
    coords = np.array([[0., 0.], [300., 400.], [600., 800.]])
    assert np.allclose(_profile_distance_km(coords), [0, .5, 1.])
    assert _is_profile_geometry(coords)
    values, depth = _smooth_depth_axis(np.array([[1.], [2.], [3.]]), np.array([0., 1., 2.]), 5)
    assert depth.size == 5 and values.shape == (5, 1)


def test_inv3d_section_plot_helpers():
    from pycsamt.agents.inv3d_agent import (
        _plot_depth_slices, _plot_resistivity_section,
        _plot_uncertainty_depth_map, _plot_uncertainty_section,
    )

    coords = np.array([[0., 0.], [500., 0.], [1000., 0.]])
    names = ["A", "B", "C"]
    depths = np.array([0., .5, 1., 2.])
    rho = np.array([[10., 20., 30., 40.], [20., 30., 40., 50.], [30., 40., 50., 60.]])
    unc = np.full_like(rho, 0.2)
    assert _plot_depth_slices(rho, coords, names, depths, 4) is not None
    assert _plot_resistivity_section(rho, names, depths, coords) is not None
    assert _plot_uncertainty_section(unc, coords, names, depths) is not None
    assert _plot_uncertainty_depth_map(unc, coords, names, depths) is not None


def test_hybrid_and_ensemble_plot_helpers():
    import pandas as pd
    from pycsamt.agents.hybrid_agent import (
        _plot_loss_curves, _plot_pinn_section, _rms_from_residuals,
    )
    from pycsamt.agents.ensemble_agent import _plot_uncertainty_section

    loss = pd.DataFrame({"epoch": [0, 1, 2], "loss": [3., 2., 1.], "data_loss": [2., 1., .5]})
    assert _plot_loss_curves(loss) is not None
    rms = _rms_from_residuals(pd.DataFrame({"residual": [3., 4.]}))
    assert isinstance(rms, tuple)
    assert _plot_pinn_section(
        np.arange(12.).reshape(3, 4), ["A", "B", "C"], 4,
        np.array([0., .5, 1., 2.]),
    ) is not None
    mean = {"A": np.arange(4.), "B": np.arange(4.) + 1}
    std = {"A": np.ones(4), "B": np.ones(4) * 2}
    assert _plot_uncertainty_section(mean, std, mean, mean, 4, np.logspace(-2, 2, 8)) is not None


def test_inv2d_array_helpers():
    from pycsamt.agents.inv2d_agent import (
        _compute_rms_2d, _inv2d_thicknesses, _nearest_station_features,
        _triangle_position_features,
    )

    assert _inv2d_thicknesses(4, np.logspace(-2, 2, 8), 1000).size == 3
    mesh = SimpleNamespace(
        nodes=np.array([[0., 0.], [1., 0.], [0., 1.]]),
        triangles=np.array([[0, 1, 2]]),
        triangle_centroids_m=np.array([[0.25, 0.25], [0.75, 0.25], [0.5, 0.75]]),
    )
    feats = [np.ones((2, 3)), np.ones((2, 3)) * 2]
    assert _nearest_station_features(feats, mesh, np.array([0., 1.])).shape[0] == 3
    assert _triangle_position_features(mesh).shape == (3, 2)
    observed = np.ones((2, 16), dtype=float)
    pred = np.ones((2, 4), dtype=float) * 2
    assert np.isnan(_compute_rms_2d(observed, pred, np.ones(3), np.logspace(-2, 2, 4)))


def test_hybrid_dimension_runners(monkeypatch):
    from pycsamt.agents.hybrid_agent import HybridInversionAgent
    import pycsamt.ai.inversion.hybrid1d as h1
    import pycsamt.ai.inversion.hybrid2d as h2
    import pycsamt.ai.inversion.hybrid3d as h3

    class Fake1D:
        def __init__(self, *args, **kwargs):
            self._ai_inv = SimpleNamespace(n_layers=2)
            self.n_sites = 2
            self._stage1 = [{"log_rho": [1, 2]}, {"log_rho": [2, 3]}]
            self._stage2 = [{"log_rho": [3, 4]}, {"log_rho": [4, 5]}]

        def fit(self, **kwargs):
            return self

        def convergence_curves(self):
            return "curve"

        def residuals(self, stage):
            return "residuals"

    class FakeND:
        def __init__(self, *args, **kwargs):
            pass

        def fit(self, **kwargs):
            return self

        def resistivity_section(self, stage=2):
            return np.ones((2, 2)) * stage

        def resistivity_volume(self, stage=2):
            return np.ones((2, 2, 2)) * stage

        def convergence_curve(self):
            return "curve"

        def residuals(self, stage):
            return "residuals"

    monkeypatch.setattr(h1, "HybridInverter1D", Fake1D)
    monkeypatch.setattr(h2, "HybridInverter2D", FakeND)
    monkeypatch.setattr(h3, "HybridInverter3D", FakeND)
    agent = HybridInversionAgent()
    for dim in (1, 2, 3):
        result = agent._run(dim, object(), object(), 2, .1, .1, .1, [])
        assert result[1] is not None and result[3] == "curve"
