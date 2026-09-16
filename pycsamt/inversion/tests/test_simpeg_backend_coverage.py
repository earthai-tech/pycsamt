# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage tests for pycsamt.inversion.backends.simpeg.SimPEGBackend.

SimPEG and discretize are not installed in this environment, and neither
is listed under pyproject.toml's ``full`` extras group, so
``test_backend_compatibility.py::test_installed_simpeg_backend_smoke``
(gated by ``importlib.util.find_spec("simpeg")``) never actually runs in
CI either -- it is a real-install smoke test only, for local use with
SimPEG manually installed.

To exercise the backend's own glue logic (mesh assembly, observation
packing, regularization/directive wiring, model conversion) without a
real SimPEG installation, this file monkeypatches the module-level
``_load_simpeg()`` loader to return light-weight fake SimPEG-shaped
modules that satisfy exactly the attribute surface ``simpeg.py`` uses.
This is the same technique the CI-realistic strategy elsewhere in this
batch relies on: none of it depends on data or packages absent from CI.
"""

from __future__ import annotations

import types

import numpy as np
import pytest

import pycsamt.inversion.backends.simpeg as simpeg_mod
from pycsamt.inversion.backends.simpeg import SimPEGBackend
from pycsamt.inversion.config import InversionConfig
from pycsamt.inversion.model import StartingModel
from pycsamt.inversion.results import InversionResult

# ---------------------------------------------------------------------------
# Fake SimPEG-shaped modules
# ---------------------------------------------------------------------------


class FakeTensorMesh:
    def __init__(self, h, origin=None):
        self.h = list(h)
        self.origin = origin
        self._shape = tuple(len(hi) for hi in self.h)
        self.nC = int(np.prod(self._shape)) if self._shape else 0

    @property
    def shape_cells(self):
        return self._shape


class FakeIdentityMap:
    def __init__(self, nP=None):
        self.nP = nP

    def __mul__(self, other):
        return FakeComposedMap(self, other)


class FakeExpMap:
    def __init__(self, mesh=None):
        self.mesh = mesh

    def __mul__(self, other):
        return FakeComposedMap(self, other)


class FakeComposedMap:
    def __init__(self, a, b):
        self.a = a
        self.b = b


class FakePointNaturalSource:
    def __init__(self, location, orientation="xy", component="apparent_resistivity"):
        self.location = location
        self.orientation = orientation
        self.component = component


class FakePlanewave:
    def __init__(self, rx_list, frequency=None):
        self.rx_list = rx_list
        self.frequency = frequency


class FakeSurvey:
    def __init__(self, source_list):
        self.source_list = source_list
        self.dobs = None


class FakeSimulation:
    def __init__(self, mesh, survey=None, sigmaMap=None, sigmaPrimary=None):
        self.mesh = mesh
        self.survey = survey
        self.sigmaMap = sigmaMap
        self.sigmaPrimary = sigmaPrimary

    def dpred(self, m):
        dobs = getattr(self.survey, "dobs", None)
        if dobs is None:
            return np.zeros(0)
        return np.asarray(dobs, dtype=float) + 0.01


class FakeData:
    def __init__(self, survey, dobs=None, standard_deviation=None):
        self.survey = survey
        self.dobs = dobs
        self.standard_deviation = standard_deviation
        if survey is not None:
            survey.dobs = dobs


class FakeDataAlwaysTypeError:
    def __init__(self, *a, **k):
        raise TypeError("simulate a SimPEG version whose Data() signature differs")


class FakeDataMisfit:
    def __init__(self, data=None, simulation=None):
        self.data = data
        self.simulation = simulation

    def __call__(self, m):
        return 0.5


class FakeRegularization:
    def __init__(self, mesh=None, mapping=None):
        self.mesh = mesh
        self.mapping = mapping
        self.alpha_s = 1.0
        self.alpha_x = 1.0
        self.alpha_y = 1.0
        self.alpha_z = 1.0
        self.reference_model = None
        self.mref = None


class FakeOptimization:
    def __init__(self, maxIter=1, maxIterLS=20, tolX=1e-5, tolF=1e-5):
        self.maxIter = maxIter
        self.maxIterLS = maxIterLS
        self.tolX = tolX
        self.tolF = tolF
        self.iter = maxIter


class FakeInvProblem:
    def __init__(self, dmis, reg, opt):
        self.dmis = dmis
        self.reg = reg
        self.opt = opt
        self.beta = 1.0


class FakeInversion:
    def __init__(self, inv_problem, directiveList=None):
        self.inv_problem = inv_problem
        self.directiveList = directiveList or []

    def run(self, m0):
        return np.asarray(m0, dtype=float) + 0.05


class FakeBetaEstimate:
    def __init__(self, beta0_ratio=1.0):
        self.beta0_ratio = beta0_ratio


class FakeTargetMisfit:
    def __init__(self, chifact=1.0):
        self.chifact = chifact


class FakeBetaSchedule:
    def __init__(self, coolingFactor=2.0, coolingRate=1):
        self.coolingFactor = coolingFactor
        self.coolingRate = coolingRate


def _make_fake_modules(
    *,
    with_planewave_xy=True,
    with_survey_cls=True,
    data_cls=None,
    with_directives=True,
    with_weighted_least_squares=True,
):
    sources_kwargs = {"Planewave": FakePlanewave}
    if with_planewave_xy:
        sources_kwargs["PlanewaveXYPrimary"] = FakePlanewave
    nsem_kwargs = dict(
        receivers=types.SimpleNamespace(PointNaturalSource=FakePointNaturalSource),
        sources=types.SimpleNamespace(**sources_kwargs),
        Simulation1DElectricField=FakeSimulation,
        Simulation3DPrimarySecondary=FakeSimulation,
        survey=types.SimpleNamespace(Data=FakeData),
    )
    if with_survey_cls:
        nsem_kwargs["Survey"] = FakeSurvey
    nsem_ns = types.SimpleNamespace(**nsem_kwargs)

    directives_kwargs = {}
    if with_directives:
        directives_kwargs = dict(
            BetaEstimate_ByEig=FakeBetaEstimate,
            TargetMisfit=FakeTargetMisfit,
            BetaSchedule=FakeBetaSchedule,
        )

    reg_kwargs = {}
    if with_weighted_least_squares:
        reg_kwargs["WeightedLeastSquares"] = FakeRegularization
    else:
        reg_kwargs["Simple"] = FakeRegularization

    return simpeg_mod._SimPEGModules(
        discretize=types.SimpleNamespace(TensorMesh=FakeTensorMesh),
        maps=types.SimpleNamespace(IdentityMap=FakeIdentityMap, ExpMap=FakeExpMap),
        data=types.SimpleNamespace(Data=data_cls or FakeData),
        data_misfit=types.SimpleNamespace(L2DataMisfit=FakeDataMisfit),
        directives=types.SimpleNamespace(**directives_kwargs),
        inverse_problem=types.SimpleNamespace(BaseInvProblem=FakeInvProblem),
        inversion=types.SimpleNamespace(BaseInversion=FakeInversion),
        optimization=types.SimpleNamespace(InexactGaussNewton=FakeOptimization),
        regularization=types.SimpleNamespace(**reg_kwargs),
        nsem=nsem_ns,
    )


@pytest.fixture
def fake_modules():
    return _make_fake_modules()


@pytest.fixture
def patch_simpeg(monkeypatch, fake_modules):
    monkeypatch.setattr(simpeg_mod, "_load_simpeg", lambda: fake_modules)
    return fake_modules


# ---------------------------------------------------------------------------
# ImportError / NotImplementedError / ValueError branches
# ---------------------------------------------------------------------------


def test_simpeg_not_installed_raises_import_error():
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0, 10.0], "rho_a": [100.0, 120.0]},
    )
    with pytest.raises(ImportError, match="SimPEG"):
        SimPEGBackend(cfg).run()


def test_load_simpeg_function_raises_import_error_directly():
    with pytest.raises(ImportError, match="SimPEG"):
        simpeg_mod._load_simpeg()


def test_unsupported_method_dimension_raises_not_implemented():
    cfg = InversionConfig(
        method="tdem",
        dimension="1d",
        backend="simpeg",
        data={"times": [1e-5, 1e-4], "values": [1e-8, 1e-9]},
    )
    with pytest.raises(NotImplementedError):
        SimPEGBackend(cfg).run()


def test_missing_mt_response_raises_value_error(patch_simpeg, tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0, 10.0]},  # no rho_a/phase
        workdir=str(tmp_path),
    )
    with pytest.raises(ValueError, match="rho_a and/or phase"):
        SimPEGBackend(cfg).run()


# ---------------------------------------------------------------------------
# 1-D sounding
# ---------------------------------------------------------------------------


def test_run_1d_sounding_success(patch_simpeg, tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={
            "freqs": [1.0, 10.0, 100.0],
            "rho_a": [100.0, 120.0, 90.0],
            "phase": [45.0, 47.0, 44.0],
        },
        max_iter=3,
        workdir=str(tmp_path / "simpeg1d"),
    )
    result = SimPEGBackend(cfg).run()

    assert isinstance(result, InversionResult)
    assert result.backend == "simpeg"
    assert result.dimension == "1d"
    assert result.status == "success"
    assert isinstance(result.model, StartingModel)
    assert result.mesh.dimension == "1d"
    assert result.metadata["engine"] == "simpeg"
    assert result.metadata["model_parameter"] == "log_sigma_cells"
    assert result.metadata["station_index"] is None
    assert result.n_iter == cfg.max_iter
    assert np.isfinite(result.rms)
    assert np.isfinite(result.objective)
    assert "mesh" in result.native and "recovered_model" in result.native


def test_run_1d_sounding_rho_only_no_phase(patch_simpeg, tmp_path):
    cfg = InversionConfig(
        method="csamt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0, 10.0], "rho_a": [80.0, 95.0]},
        max_iter=2,
        workdir=str(tmp_path / "simpeg1d_rho_only"),
    )
    result = SimPEGBackend(cfg).run()
    assert result.status == "success"


def test_run_1d_sounding_custom_backend_options(patch_simpeg, tmp_path):
    cfg = InversionConfig(
        method="amt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0, 10.0], "rho_a": [80.0, 95.0], "phase": [40.0, 42.0]},
        regularization="damped",
        backend_options={
            "estimate_beta": False,
            "beta0": 2.0,
            "target_chifact": 1.5,
            "cooling_factor": 3.0,
            "cooling_rate": 2,
            "max_iter_ls": 5,
        },
        max_iter=2,
        workdir=str(tmp_path / "simpeg1d_opts"),
    )
    result = SimPEGBackend(cfg).run()
    assert result.metadata["beta0"] == 2.0


# ---------------------------------------------------------------------------
# 2-D stitched profile
# ---------------------------------------------------------------------------


def test_run_2d_profile_all_stations_succeed(patch_simpeg, tmp_path):
    cfg = InversionConfig(
        method="amt",
        dimension="2d",
        backend="simpeg",
        data={
            "freqs": [10.0, 100.0],
            "rho_a": [[80.0, 100.0], [90.0, 110.0]],
            "phase": [[42.0, 45.0], [43.0, 46.0]],
            "station_x": [0.0, 250.0],
            "station_names": ["A01", "A02"],
        },
        max_iter=2,
        workdir=str(tmp_path / "simpeg2d"),
    )
    result = SimPEGBackend(cfg).run()

    assert result.dimension == "2d"
    assert result.status == "success"
    assert result.warnings == []
    assert result.model["rho_2d"].shape[1] == 2
    assert result.metadata["profile_mode"] == "stitched_station_1d"
    assert len(result.native) == 2


def test_run_2d_profile_one_station_fails(patch_simpeg, monkeypatch, tmp_path):
    cfg = InversionConfig(
        method="amt",
        dimension="2d",
        backend="simpeg",
        data={
            "freqs": [10.0, 100.0],
            "rho_a": [[80.0, 100.0], [90.0, 110.0]],
            "phase": [[42.0, 45.0], [43.0, 46.0]],
            "station_x": [0.0, 250.0],
            "station_names": ["A01", "A02"],
        },
        max_iter=2,
        workdir=str(tmp_path / "simpeg2d_partial"),
    )
    original = SimPEGBackend._run_sounding

    def _flaky(self, em_data, modules, *, station_index):
        if station_index == 1:
            raise RuntimeError("synthetic station failure")
        return original(self, em_data, modules, station_index=station_index)

    monkeypatch.setattr(SimPEGBackend, "_run_sounding", _flaky)
    result = SimPEGBackend(cfg).run()

    assert result.status == "needs_review"
    assert len(result.warnings) == 1
    assert "A02" in result.warnings[0]
    assert result.model["station_names"] == ["A01"]


def test_run_2d_profile_all_stations_fail_raises(patch_simpeg, monkeypatch, tmp_path):
    cfg = InversionConfig(
        method="amt",
        dimension="2d",
        backend="simpeg",
        data={
            "freqs": [10.0],
            "rho_a": [[80.0], [90.0]],
            "station_x": [0.0, 250.0],
        },
        max_iter=1,
        workdir=str(tmp_path / "simpeg2d_all_fail"),
    )

    def _always_fails(self, em_data, modules, *, station_index):
        raise RuntimeError("synthetic total failure")

    monkeypatch.setattr(SimPEGBackend, "_run_sounding", _always_fails)
    with pytest.raises(RuntimeError, match="all SimPEG station inversions failed"):
        SimPEGBackend(cfg).run()


# ---------------------------------------------------------------------------
# 3-D primary-secondary
# ---------------------------------------------------------------------------


def test_run_3d_success(patch_simpeg, tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="simpeg",
        data={
            "freqs": [1.0],
            "rho_a": [[100.0], [120.0]],
            "phase": [[45.0], [47.0]],
            "station_x": [0.0, 500.0],
            "station_names": ["M01", "M02"],
            "metadata": {"station_y": [0.0, 250.0]},
        },
        backend_options={
            "sigma_primary": 0.01,
            "beta0": 1.0,
            "target_chifact": 1.0,
            "nx": 4,
            "ny": 4,
            "nz": 4,
        },
        max_iter=2,
        workdir=str(tmp_path / "simpeg3d"),
    )
    result = SimPEGBackend(cfg).run()

    assert result.dimension == "3d"
    assert result.status == "success"
    assert result.model["rho_3d"].shape == (4, 4, 4)
    assert result.mesh.metadata["mesh_shape"] == (4, 4, 4)
    assert result.metadata["sigma_primary"] == 0.01
    assert result.metadata["simulation"] == "Simulation3DPrimarySecondary"


def test_run_3d_without_planewave_xy_primary_falls_back_to_planewave(
    monkeypatch, tmp_path
):
    modules = _make_fake_modules(with_planewave_xy=False)
    monkeypatch.setattr(simpeg_mod, "_load_simpeg", lambda: modules)
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="simpeg",
        data={
            "freqs": [1.0],
            "rho_a": [[100.0], [120.0]],
            "station_x": [0.0, 500.0],
        },
        backend_options={"nx": 2, "ny": 2, "nz": 2},
        max_iter=1,
        workdir=str(tmp_path / "simpeg3d_no_xy"),
    )
    result = SimPEGBackend(cfg).run()
    assert result.status == "success"


def test_run_3d_survey_fallback_import_path(monkeypatch, tmp_path):
    """When the installed nsem module has no ``Survey`` attribute, the
    backend falls back to ``from simpeg import survey as survey_mod``."""
    modules = _make_fake_modules(with_survey_cls=False)
    fake_survey_mod = types.SimpleNamespace(Survey=FakeSurvey)
    fake_simpeg_pkg = types.ModuleType("simpeg")
    fake_simpeg_pkg.survey = fake_survey_mod
    monkeypatch.setitem(
        __import__("sys").modules, "simpeg", fake_simpeg_pkg
    )
    monkeypatch.setitem(
        __import__("sys").modules, "simpeg.survey", fake_survey_mod
    )
    monkeypatch.setattr(simpeg_mod, "_load_simpeg", lambda: modules)

    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0, 10.0], "rho_a": [100.0, 120.0]},
        max_iter=1,
        workdir=str(tmp_path / "simpeg_survey_fallback"),
    )
    result = SimPEGBackend(cfg).run()
    assert result.status == "success"


# ---------------------------------------------------------------------------
# _build_simpeg_data fallback
# ---------------------------------------------------------------------------


def test_build_simpeg_data_falls_back_on_type_error():
    modules = _make_fake_modules(data_cls=FakeDataAlwaysTypeError)
    survey = FakeSurvey([])
    observed = np.array([1.0, 2.0])
    errors = np.array([0.1, 0.1])
    result = simpeg_mod._build_simpeg_data(survey, observed, errors, modules)
    assert isinstance(result, FakeData)
    assert result.dobs is observed


def test_build_simpeg_data_normal_path():
    modules = _make_fake_modules()
    survey = FakeSurvey([])
    observed = np.array([1.0, 2.0])
    errors = np.array([0.1, 0.1])
    result = simpeg_mod._build_simpeg_data(survey, observed, errors, modules)
    assert isinstance(result, FakeData)


# ---------------------------------------------------------------------------
# _build_regularization
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kind", ["none", "damped", "smooth"])
def test_build_regularization_kinds(kind):
    modules = _make_fake_modules()
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
        regularization=kind,
    )
    mesh = FakeTensorMesh([np.array([1.0, 2.0])])
    mapping = FakeIdentityMap(nP=2)
    reg = simpeg_mod._build_regularization(mesh, mapping, modules, cfg)

    if kind == "none":
        assert reg.alpha_s == 0.0 and reg.alpha_x == 0.0
        assert reg.alpha_y == 0.0 and reg.alpha_z == 0.0
    elif kind == "damped":
        assert reg.alpha_s == 1.0
        assert reg.alpha_x == 0.0 and reg.alpha_y == 0.0 and reg.alpha_z == 0.0
    else:
        assert reg.alpha_s == 1.0
        assert reg.alpha_x == 1.0
        assert reg.alpha_z == 1.0


def test_build_regularization_uses_simple_when_no_weighted_least_squares():
    modules = _make_fake_modules(with_weighted_least_squares=False)
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
    )
    mesh = FakeTensorMesh([np.array([1.0, 2.0])])
    mapping = FakeIdentityMap(nP=2)
    reg = simpeg_mod._build_regularization(mesh, mapping, modules, cfg)
    assert isinstance(reg, FakeRegularization)


def test_build_regularization_sets_reference_model_when_present():
    modules = _make_fake_modules()
    mesh = FakeTensorMesh([np.array([1.0, 2.0])])
    mapping = FakeIdentityMap(nP=2)
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
        reference_model=StartingModel([100.0, 100.0], [50.0]),
    )
    reg = simpeg_mod._build_regularization(mesh, mapping, modules, cfg)
    assert reg.reference_model is not None
    assert reg.mref is not None


# ---------------------------------------------------------------------------
# _simpeg_reference_model
# ---------------------------------------------------------------------------


def test_simpeg_reference_model_none_when_unset():
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
    )
    mesh = FakeTensorMesh([np.array([1.0, 2.0])])
    assert simpeg_mod._simpeg_reference_model(cfg, mesh) is None


def test_simpeg_reference_model_positive_resistivity_converted_to_log_sigma():
    mesh = FakeTensorMesh([np.array([1.0, 2.0])])
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
        reference_model=StartingModel([50.0, 200.0], [10.0]),
    )
    out = simpeg_mod._simpeg_reference_model(cfg, mesh)
    np.testing.assert_allclose(out, np.log(1.0 / np.array([50.0, 200.0])))


def test_simpeg_reference_model_size_mismatch_returns_none():
    mesh = FakeTensorMesh([np.array([1.0, 2.0, 3.0])])  # nC = 3
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
        reference_model=StartingModel([50.0, 200.0], [10.0]),  # size 2
    )
    assert simpeg_mod._simpeg_reference_model(cfg, mesh) is None


def test_simpeg_reference_model_backend_options_override():
    mesh = FakeTensorMesh([np.array([1.0, 2.0])])
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
        backend_options={"reference_model": np.array([-1.0, -2.0])},
    )
    out = simpeg_mod._simpeg_reference_model(cfg, mesh)
    np.testing.assert_allclose(out, np.array([-1.0, -2.0]))


def test_simpeg_reference_model_unconvertible_returns_none():
    mesh = FakeTensorMesh([np.array([1.0, 2.0])])
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
        backend_options={"reference_model": object()},
    )
    assert simpeg_mod._simpeg_reference_model(cfg, mesh) is None


# ---------------------------------------------------------------------------
# _build_directives
# ---------------------------------------------------------------------------


def test_build_directives_defaults_returns_three():
    modules = _make_fake_modules()
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
    )
    directives = simpeg_mod._build_directives(modules, cfg)
    assert len(directives) == 3


def test_build_directives_estimate_beta_false_skips_beta_directive():
    modules = _make_fake_modules()
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
        backend_options={"estimate_beta": False},
    )
    directives = simpeg_mod._build_directives(modules, cfg)
    assert len(directives) == 2
    assert not any(isinstance(d, FakeBetaEstimate) for d in directives)


def test_build_directives_empty_namespace_returns_empty_list():
    modules = _make_fake_modules(with_directives=False)
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0]},
    )
    assert simpeg_mod._build_directives(modules, cfg) == []


# ---------------------------------------------------------------------------
# _mesh_shape
# ---------------------------------------------------------------------------


def test_mesh_shape_from_shape_cells():
    mesh = FakeTensorMesh([np.ones(2), np.ones(3), np.ones(4)])
    assert simpeg_mod._mesh_shape(mesh) == (2, 3, 4)


def test_mesh_shape_from_vnc_attribute():
    mesh = types.SimpleNamespace(vnC=(2, 3, 4))
    assert simpeg_mod._mesh_shape(mesh) == (2, 3, 4)


def test_mesh_shape_wrong_ndim_raises_value_error():
    mesh = types.SimpleNamespace(vnC=(2, 3))
    with pytest.raises(ValueError, match="3-D mesh"):
        simpeg_mod._mesh_shape(mesh)


def test_mesh_shape_missing_attrs_raises_attribute_error():
    mesh = types.SimpleNamespace()
    with pytest.raises(AttributeError):
        simpeg_mod._mesh_shape(mesh)


# ---------------------------------------------------------------------------
# starting/recovered model conversions
# ---------------------------------------------------------------------------


def test_starting_sigma_model_matches_layer_by_depth():
    start = StartingModel([10.0, 100.0], [50.0])
    z_centers = np.array([10.0, 60.0])
    out = simpeg_mod._starting_sigma_model(start, z_centers)
    np.testing.assert_allclose(out[0], np.log(1.0 / 10.0))
    np.testing.assert_allclose(out[1], np.log(1.0 / 100.0))


def test_starting_3d_log_sigma_shape_and_fortran_order():
    start = StartingModel([10.0, 100.0], [50.0])
    centers = {
        "x": np.arange(2.0),
        "y": np.arange(3.0),
        "z_depth": np.array([10.0, 60.0]),
        "z": np.array([-10.0, -60.0]),
    }
    out = simpeg_mod._starting_3d_log_sigma(start, centers)
    assert out.shape == (2 * 3 * 2,)


def test_rho_3d_from_log_sigma_roundtrip():
    mesh = FakeTensorMesh([np.ones(2), np.ones(2), np.ones(2)])
    log_sigma = np.log(1.0 / np.array([10.0] * 8))
    rho = simpeg_mod._rho_3d_from_log_sigma(log_sigma, mesh)
    assert rho.shape == (2, 2, 2)
    np.testing.assert_allclose(rho, 10.0, rtol=1e-6)


def test_model_from_sigma_cells_normal_and_empty_bin_fallback():
    z_centers = np.array([0.0, 100.0])
    log_sigma = np.array([np.log(1.0 / 10.0), np.log(1.0 / 1000.0)])
    recovered = simpeg_mod._model_from_sigma_cells(log_sigma, z_centers, 5)
    assert recovered.n_layers == 5
    assert recovered.thicknesses.size == 4
    assert np.all(recovered.resistivities > 0)


def test_model_from_sigma_cells_clamps_n_layers_below_two():
    z_centers = np.array([0.0, 50.0, 100.0])
    log_sigma = np.log(1.0 / np.array([10.0, 20.0, 30.0]))
    recovered = simpeg_mod._model_from_sigma_cells(log_sigma, z_centers, 1)
    assert recovered.n_layers == 2


def test_layer_centers_normal():
    thicknesses = np.array([10.0, 20.0])
    centers = simpeg_mod._layer_centers(thicknesses)
    assert centers.size == 3
    assert centers[0] == 5.0


def test_layer_centers_empty_thicknesses():
    centers = simpeg_mod._layer_centers(np.array([]))
    assert centers.size == 1
    assert centers[0] == 0.5


# ---------------------------------------------------------------------------
# station helpers
# ---------------------------------------------------------------------------


def test_station_helpers_defaults_and_overrides():
    from pycsamt.inversion.data import EMData

    em_data_defaults = EMData(
        method="mt",
        frequencies=[1.0, 10.0],
        rho_a=[[80.0, 90.0], [70.0, 60.0]],
    )
    names = simpeg_mod._station_names(em_data_defaults, 2)
    assert names == ["S000", "S001"]
    xs = simpeg_mod._station_x(em_data_defaults, 2)
    np.testing.assert_allclose(xs, [0.0, 1.0])
    ys = simpeg_mod._station_y(em_data_defaults, 2)
    np.testing.assert_allclose(ys, [0.0, 0.0])

    em_data_named = EMData(
        method="mt",
        frequencies=[1.0, 10.0],
        rho_a=[[80.0, 90.0], [70.0, 60.0]],
        station_names=["X1", "X2"],
        station_x=[5.0, 15.0],
        metadata={"station_y": [1.0, 2.0]},
    )
    assert simpeg_mod._station_names(em_data_named, 2) == ["X1", "X2"]
    np.testing.assert_allclose(
        simpeg_mod._station_x(em_data_named, 2), [5.0, 15.0]
    )
    np.testing.assert_allclose(
        simpeg_mod._station_y(em_data_named, 2), [1.0, 2.0]
    )

    locs = simpeg_mod._station_locations(em_data_named)
    assert locs.shape == (2, 3)
    np.testing.assert_allclose(locs[:, 2], [0.0, 0.0])


@pytest.mark.parametrize("meta_key", ["station_y", "y", "northing"])
def test_station_y_reads_all_recognized_metadata_keys(meta_key):
    from pycsamt.inversion.data import EMData

    em_data = EMData(
        method="mt",
        frequencies=[1.0],
        rho_a=[[80.0], [90.0]],
        metadata={meta_key: [3.0, 4.0]},
    )
    np.testing.assert_allclose(simpeg_mod._station_y(em_data, 2), [3.0, 4.0])


def test_row_1d_and_2d_branches():
    arr_1d = np.array([1.0, 2.0, 3.0])
    np.testing.assert_allclose(simpeg_mod._row(arr_1d, 1), arr_1d)
    arr_2d = np.array([[1.0, 2.0], [3.0, 4.0]])
    np.testing.assert_allclose(simpeg_mod._row(arr_2d, 1), [3.0, 4.0])


def test_station_y_metadata_size_mismatch_falls_through_to_zeros():
    from pycsamt.inversion.data import EMData

    em_data = EMData(
        method="mt",
        frequencies=[1.0],
        rho_a=[[80.0], [90.0]],
        metadata={"station_y": [1.0]},  # size 1, but n_st=2 requested
    )
    np.testing.assert_allclose(simpeg_mod._station_y(em_data, 2), [0.0, 0.0])


def test_run_1d_sounding_phase_only_skips_rho_a_receiver(patch_simpeg, tmp_path):
    """Exercises the ``em_data.rho_a is None`` branch of
    ``_build_nsem_survey`` and ``_pack_nsem_observations`` (phase-only
    natural-source data, no apparent resistivity)."""
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0, 10.0], "phase": [42.0, 44.0]},
        max_iter=1,
        workdir=str(tmp_path / "simpeg_phase_only"),
    )
    result = SimPEGBackend(cfg).run()
    assert result.status == "success"


def test_run_2d_profile_missing_mesh_falls_back_to_arange_z_centers(
    patch_simpeg, monkeypatch, tmp_path
):
    """When every stitched-station result carries ``mesh=None``, the
    profile assembly falls back to ``np.arange`` depth centers."""
    from pycsamt.inversion.results import InversionResult as _Result

    cfg = InversionConfig(
        method="amt",
        dimension="2d",
        backend="simpeg",
        data={
            "freqs": [10.0, 100.0],
            "rho_a": [[80.0, 100.0], [90.0, 110.0]],
            "station_x": [0.0, 250.0],
            "station_names": ["A01", "A02"],
        },
        max_iter=1,
        workdir=str(tmp_path / "simpeg2d_no_mesh"),
    )

    def _no_mesh_sounding(self, em_data, modules, *, station_index):
        return _Result(
            method=cfg.method,
            dimension="1d",
            backend="simpeg",
            status="success",
            model=StartingModel([100.0, 100.0], [50.0]),
            mesh=None,
            rms=1.0,
            objective=0.1,
            n_iter=1,
        )

    monkeypatch.setattr(SimPEGBackend, "_run_sounding", _no_mesh_sounding)
    result = SimPEGBackend(cfg).run()
    assert result.status == "success"
    np.testing.assert_allclose(result.mesh.z_centers, [0.0, 1.0])


def test_station_data_builds_single_station_emdata():
    from pycsamt.inversion.data import EMData

    em_data = EMData(
        method="mt",
        frequencies=[1.0, 10.0],
        rho_a=[[80.0, 90.0], [70.0, 60.0]],
        phase=[[40.0, 41.0], [42.0, 43.0]],
        errors=[[1.0, 1.0], [2.0, 2.0]],
        station_names=["A", "B"],
        station_x=[0.0, 100.0],
    )
    single = simpeg_mod._station_data(em_data, 1)
    assert single.station_names == ["B"]
    np.testing.assert_allclose(single.rho_a, [70.0, 60.0])
    np.testing.assert_allclose(single.phase, [42.0, 43.0])
    np.testing.assert_allclose(single.station_x, [100.0])
