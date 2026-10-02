# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Coverage tests for pycsamt.inversion.backends.pygimli.PyGIMLiBackend.

pyGIMLi is not installed in this environment (``python -c "import pygimli"``
raises ``ModuleNotFoundError``) and is not listed under pyproject.toml's
``full`` extras group, so
``test_backend_compatibility.py::test_installed_pygimli_backend_smoke``
(gated by ``importlib.util.find_spec("pygimli")``) never runs in CI either --
it is a real-install smoke test only. ``test_inversion_api.py`` already
exercises several of this module's standalone helper functions
(``_pack_mt_observations/_errors``, the MT/TDEM operator-fallback
constructors, ``_run_inversion``'s relative-error fallback) directly and the
real ``ImportError`` path; this file fills the remaining, larger gap: the
end-to-end ``PyGIMLiBackend.run()`` flow (1-D MT/TDEM soundings, stitched
2-D profiles) plus the smaller untested helpers, all via the same
monkeypatched-fake-module technique used for ``simpeg.py`` in this batch --
none of it depends on data or packages absent from CI.
"""

from __future__ import annotations

import numpy as np
import pytest

import pycsamt.inversion.backends.pygimli as pygimli_mod
from pycsamt.inversion.backends.pygimli import PyGIMLiBackend
from pycsamt.inversion.config import InversionConfig
from pycsamt.inversion.data import EMData
from pycsamt.inversion.model import StartingModel
from pycsamt.inversion.results import InversionResult

# ---------------------------------------------------------------------------
# Fake pyGIMLi-shaped modules
# ---------------------------------------------------------------------------

_RESPONSE_SIZE = {"n": 2}


class FakeMTOperator:
    """Stand-in for MT1dSmoothModelling / MT1DSmoothModelling."""

    def __init__(self, T=None, thk=None, periods=None, verbose=False, **_extra):
        self.T = T if T is not None else periods
        self.thk = thk

    def response(self, model):
        return np.linspace(1.0, 2.0, _RESPONSE_SIZE["n"])


class FakeMTBlockOperator:
    """Stand-in for MT1dBlockModelling: only accepts nLayers-style kwargs."""

    def __init__(self, T=None, nLayers=None, verbose=False, **_extra):
        if nLayers is None:
            raise TypeError("nLayers required")
        self.T = T
        self.nLayers = nLayers

    def response(self, model):
        return np.linspace(1.0, 2.0, _RESPONSE_SIZE["n"])


class FakeTDEMOperator:
    def __init__(
        self, thk=None, times=None, txArea=None, rxArea=None, verbose=False, **_extra
    ):
        self.thk = thk
        self.times = times

    def response(self, model):
        return np.full(len(self.times), 1e-8)


class FakeInversion:
    def __init__(self, fop=None, verbose=False):
        self.fop = fop
        self.verbose = verbose
        self.chi2 = 1.2
        self.iter = 4

    def run(self, observed, **kwargs):
        start = np.asarray(kwargs["startModel"], dtype=float)
        return start * 1.01 + 0.01


class FakeTrans:
    class TransLog:
        pass


def _make_fake_em(*, block=False):
    if block:
        return type(
            "FakeEM",
            (),
            {"MT1dBlockModelling": FakeMTBlockOperator, "TDEMSmoothModelling": FakeTDEMOperator},
        )
    return type(
        "FakeEM",
        (),
        {"MT1dSmoothModelling": FakeMTOperator, "TDEMSmoothModelling": FakeTDEMOperator},
    )


def _make_fake_pg(*, with_trans=True):
    attrs = {"Inversion": FakeInversion}
    if with_trans:
        attrs["trans"] = FakeTrans
    return type("FakePG", (), attrs)


@pytest.fixture
def fake_modules():
    return pygimli_mod._PyGIMLiModules(pg=_make_fake_pg(), em=_make_fake_em())


@pytest.fixture
def patch_pygimli(monkeypatch, fake_modules):
    monkeypatch.setattr(pygimli_mod, "_load_pygimli", lambda: fake_modules)
    return fake_modules


# ---------------------------------------------------------------------------
# check_supported / ValueError guards
# ---------------------------------------------------------------------------


def test_unsupported_dimension_raises_not_implemented():
    cfg = InversionConfig(
        method="mt",
        dimension="3d",
        backend="pygimli",
        data={"freqs": [1.0], "rho_a": [100.0]},
    )
    with pytest.raises(NotImplementedError):
        PyGIMLiBackend(cfg).run()


def test_missing_tdem_response_raises_value_error(patch_pygimli, tmp_path):
    cfg = InversionConfig(
        method="tdem",
        dimension="1d",
        backend="pygimli",
        data={"times": [1e-5, 1e-4]},  # no values
        workdir=str(tmp_path),
    )
    with pytest.raises(ValueError, match="times plus values"):
        PyGIMLiBackend(cfg).run()


def test_missing_mt_response_raises_value_error(patch_pygimli, tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="pygimli",
        data={"freqs": [1.0, 10.0]},  # no rho_a/phase
        workdir=str(tmp_path),
    )
    with pytest.raises(ValueError, match="rho_a and/or phase"):
        PyGIMLiBackend(cfg).run()


# ---------------------------------------------------------------------------
# 1-D MT sounding -- full run()
# ---------------------------------------------------------------------------


def test_run_1d_mt_sounding_success(patch_pygimli, tmp_path):
    _RESPONSE_SIZE["n"] = 4
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="pygimli",
        data={
            "freqs": [1.0, 10.0],
            "rho_a": [100.0, 120.0],
            "phase": [45.0, 47.0],
        },
        max_iter=3,
        workdir=str(tmp_path / "pg1d"),
    )
    result = PyGIMLiBackend(cfg).run()

    assert isinstance(result, InversionResult)
    assert result.backend == "pygimli"
    assert result.dimension == "1d"
    assert result.status == "success"
    assert isinstance(result.model, StartingModel)
    assert result.mesh.dimension == "1d"
    assert result.metadata["engine"] == "pygimli"
    assert result.metadata["mode"] == "mt"
    assert result.metadata["station_index"] is None
    assert result.n_iter == 4
    assert result.objective == pytest.approx(1.2)
    assert np.isfinite(result.rms)
    assert result.native["fop"] is not None


def test_run_1d_mt_sounding_no_trans_attribute(monkeypatch, tmp_path):
    """Exercises the ``hasattr(pg, "trans")`` False branch."""
    _RESPONSE_SIZE["n"] = 4
    modules = pygimli_mod._PyGIMLiModules(
        pg=_make_fake_pg(with_trans=False), em=_make_fake_em()
    )
    monkeypatch.setattr(pygimli_mod, "_load_pygimli", lambda: modules)
    cfg = InversionConfig(
        method="amt",
        dimension="1d",
        backend="pygimli",
        data={"freqs": [1.0, 10.0], "rho_a": [80.0, 95.0], "phase": [40.0, 42.0]},
        max_iter=2,
        workdir=str(tmp_path / "pg1d_notrans"),
    )
    result = PyGIMLiBackend(cfg).run()
    assert result.status == "success"


def test_run_1d_mt_sounding_block_operator(monkeypatch, tmp_path):
    """Block-model operator path -> thickness_resistivity parameterization."""
    _RESPONSE_SIZE["n"] = 2
    modules = pygimli_mod._PyGIMLiModules(pg=_make_fake_pg(), em=_make_fake_em(block=True))
    monkeypatch.setattr(pygimli_mod, "_load_pygimli", lambda: modules)
    cfg = InversionConfig(
        method="csamt",
        dimension="1d",
        backend="pygimli",
        data={"freqs": [1.0, 10.0], "rho_a": [80.0, 95.0]},
        backend_options={"mt_operator": "MT1dBlockModelling"},
        max_iter=2,
        workdir=str(tmp_path / "pg1d_block"),
    )
    result = PyGIMLiBackend(cfg).run()
    assert result.status == "success"
    assert result.metadata["parameterization"] == "thickness_resistivity"


# ---------------------------------------------------------------------------
# 1-D TDEM sounding -- full run()
# ---------------------------------------------------------------------------


def test_run_1d_tdem_sounding_success(patch_pygimli, tmp_path):
    cfg = InversionConfig(
        method="tdem",
        dimension="1d",
        backend="pygimli",
        data={"times": [1e-5, 3e-5, 1e-4], "values": [1e-8, 5e-9, 1e-9]},
        backend_options={"tx_area": 7850.0, "rx_area": 100.0},
        max_iter=3,
        workdir=str(tmp_path / "pgtdem"),
    )
    result = PyGIMLiBackend(cfg).run()
    assert result.status == "success"
    assert result.metadata["mode"] == "tdem"
    assert result.metadata["tx_area"] == 7850.0
    assert result.metadata["rx_area"] == 100.0


# ---------------------------------------------------------------------------
# 2-D stitched profile
# ---------------------------------------------------------------------------


def test_run_2d_profile_all_stations_succeed(patch_pygimli, tmp_path):
    _RESPONSE_SIZE["n"] = 4
    cfg = InversionConfig(
        method="amt",
        dimension="2d",
        backend="pygimli",
        data={
            "freqs": [10.0, 100.0],
            "rho_a": [[80.0, 100.0], [90.0, 110.0]],
            "phase": [[42.0, 45.0], [43.0, 46.0]],
            "station_x": [0.0, 250.0],
            "station_names": ["A01", "A02"],
        },
        max_iter=2,
        workdir=str(tmp_path / "pg2d"),
    )
    result = PyGIMLiBackend(cfg).run()

    assert result.dimension == "2d"
    assert result.status == "success"
    assert result.warnings == []
    assert result.model["rho_2d"].shape[1] == 2
    assert result.metadata["profile_mode"] == "stitched_station_1d"
    assert len(result.native) == 2


def test_run_2d_profile_one_station_fails(patch_pygimli, monkeypatch, tmp_path):
    _RESPONSE_SIZE["n"] = 4
    cfg = InversionConfig(
        method="amt",
        dimension="2d",
        backend="pygimli",
        data={
            "freqs": [10.0, 100.0],
            "rho_a": [[80.0, 100.0], [90.0, 110.0]],
            "phase": [[42.0, 45.0], [43.0, 46.0]],
            "station_x": [0.0, 250.0],
            "station_names": ["A01", "A02"],
        },
        max_iter=2,
        workdir=str(tmp_path / "pg2d_partial"),
    )
    original = PyGIMLiBackend._run_sounding

    def _flaky(self, em_data, modules, *, station_index):
        if station_index == 1:
            raise RuntimeError("synthetic station failure")
        return original(self, em_data, modules, station_index=station_index)

    monkeypatch.setattr(PyGIMLiBackend, "_run_sounding", _flaky)
    result = PyGIMLiBackend(cfg).run()

    assert result.status == "needs_review"
    assert len(result.warnings) == 1
    assert "A02" in result.warnings[0]
    assert result.model["station_names"] == ["A01"]


def test_run_2d_profile_all_stations_fail_raises(patch_pygimli, monkeypatch, tmp_path):
    cfg = InversionConfig(
        method="amt",
        dimension="2d",
        backend="pygimli",
        data={
            "freqs": [10.0],
            "rho_a": [[80.0], [90.0]],
            "station_x": [0.0, 250.0],
        },
        max_iter=1,
        workdir=str(tmp_path / "pg2d_all_fail"),
    )

    def _always_fails(self, em_data, modules, *, station_index):
        raise RuntimeError("synthetic total failure")

    monkeypatch.setattr(PyGIMLiBackend, "_run_sounding", _always_fails)
    with pytest.raises(RuntimeError, match="all pyGIMLi station inversions failed"):
        PyGIMLiBackend(cfg).run()


def test_run_2d_profile_missing_mesh_falls_back_to_arange_z_centers(
    patch_pygimli, monkeypatch, tmp_path
):
    cfg = InversionConfig(
        method="amt",
        dimension="2d",
        backend="pygimli",
        data={
            "freqs": [10.0, 100.0],
            "rho_a": [[80.0, 100.0], [90.0, 110.0]],
            "station_x": [0.0, 250.0],
            "station_names": ["A01", "A02"],
        },
        max_iter=1,
        workdir=str(tmp_path / "pg2d_no_mesh"),
    )

    def _no_mesh_sounding(self, em_data, modules, *, station_index):
        return InversionResult(
            method=cfg.method,
            dimension="1d",
            backend="pygimli",
            status="success",
            model=StartingModel([100.0, 100.0], [50.0]),
            mesh=None,
            rms=1.0,
            objective=0.1,
            n_iter=1,
        )

    monkeypatch.setattr(PyGIMLiBackend, "_run_sounding", _no_mesh_sounding)
    result = PyGIMLiBackend(cfg).run()
    assert result.status == "success"
    np.testing.assert_allclose(result.mesh.z_centers, [0.0, 1.0])


# ---------------------------------------------------------------------------
# _resolve_operator
# ---------------------------------------------------------------------------


def test_resolve_operator_found():
    em = _make_fake_em()
    cls, name = pygimli_mod._resolve_operator(
        em, ("Missing", "MT1dSmoothModelling"), "MT/AMT/CSAMT 1-D"
    )
    assert cls is FakeMTOperator
    assert name == "MT1dSmoothModelling"


def test_resolve_operator_not_found_raises_not_implemented():
    em = type("EmptyEM", (), {})
    with pytest.raises(NotImplementedError, match="does not expose"):
        pygimli_mod._resolve_operator(em, ("Nope1", "Nope2"), "MT/AMT/CSAMT 1-D")


# ---------------------------------------------------------------------------
# _construct_operator
# ---------------------------------------------------------------------------


def test_construct_operator_all_attempts_fail_raises_type_error():
    class AlwaysFails:
        def __init__(self, **kwargs):
            raise TypeError("nope")

    with pytest.raises(TypeError, match="could not construct pyGIMLi operator"):
        pygimli_mod._construct_operator(
            AlwaysFails, [{"a": 1}, {"b": 2}], "AlwaysFails"
        )


# ---------------------------------------------------------------------------
# _make_inversion
# ---------------------------------------------------------------------------


def test_make_inversion_kwarg_signature():
    pg = _make_fake_pg()
    inv = pygimli_mod._make_inversion(pg, fop="FOP", verbose=True)
    assert isinstance(inv, FakeInversion)
    assert inv.fop == "FOP"


def test_make_inversion_positional_fallback():
    class PositionalOnlyInversion:
        def __init__(self, fop, verbose=False):
            self.fop = fop
            self.verbose = verbose

    class PgPositional:
        Inversion = staticmethod(
            lambda *a, **k: (
                PositionalOnlyInversion(*a)
                if not k or "fop" not in k
                else (_ for _ in ()).throw(TypeError("no kwarg fop"))
            )
        )

    inv = pygimli_mod._make_inversion(PgPositional, fop="FOP", verbose=False)
    assert isinstance(inv, PositionalOnlyInversion)
    assert inv.fop == "FOP"


# ---------------------------------------------------------------------------
# _run_inversion -- minimal-kwargs final fallback
# ---------------------------------------------------------------------------


def test_run_inversion_minimal_kwargs_fallback():
    class PickyInversion:
        def run(self, observed, **kwargs):
            if set(kwargs) != {"startModel", "lam", "maxIter", "verbose"}:
                raise TypeError("only accepts the minimal signature")
            return np.asarray(kwargs["startModel"], dtype=float) + 5.0

    recovered = pygimli_mod._run_inversion(
        PickyInversion(),
        np.array([1.0, 2.0]),
        startModel=np.array([1.0, 2.0]),
        errorVals=np.array([0.1, 0.1]),
        lam=10.0,
        maxIter=2,
        verbose=False,
    )
    np.testing.assert_allclose(recovered, [6.0, 7.0])


# ---------------------------------------------------------------------------
# _recover_mt_model
# ---------------------------------------------------------------------------


def test_recover_mt_model_exact_n_layers():
    start = StartingModel([100.0, 200.0], [50.0])
    raw = np.array([10.0, 20.0])
    out = pygimli_mod._recover_mt_model(raw, start)
    np.testing.assert_allclose(out.resistivities, raw)
    np.testing.assert_allclose(out.thicknesses, start.thicknesses)


def test_recover_mt_model_block_layout():
    start = StartingModel([100.0, 200.0, 300.0], [50.0, 150.0])
    raw = np.array([10.0, 20.0, 1.0, 2.0, 3.0])  # thk(2) + rho(3)
    out = pygimli_mod._recover_mt_model(raw, start)
    np.testing.assert_allclose(out.thicknesses, [10.0, 20.0])
    np.testing.assert_allclose(out.resistivities, [1.0, 2.0, 3.0])


def test_recover_mt_model_fallback_layout():
    start = StartingModel([100.0, 200.0], [50.0])
    raw = np.array([1.0, 2.0, 3.0, 4.0, 5.0])  # unexpected size
    out = pygimli_mod._recover_mt_model(raw, start)
    np.testing.assert_allclose(out.resistivities, raw[-2:])
    np.testing.assert_allclose(out.thicknesses, start.thicknesses)


# ---------------------------------------------------------------------------
# _result_from_pygimli -- chi2/iter attribute variants
# ---------------------------------------------------------------------------


def test_result_from_pygimli_chi2_as_plain_attribute():
    class Inv:
        chi2 = 0.7
        iter = 3

    cfg = InversionConfig(
        method="mt", dimension="1d", backend="pygimli",
        data={"freqs": [1.0], "rho_a": [100.0]},
    )
    em_data = EMData.coerce({"freqs": [1.0], "rho_a": [100.0]})
    recovered = StartingModel([100.0, 200.0], [50.0])
    result = pygimli_mod._result_from_pygimli(
        cfg, em_data, recovered, np.array([100.0]), 0.1,
        Inv(), fop="FOP", station_index=0, extra={"mode": "mt"},
    )
    assert result.objective == pytest.approx(0.7)
    assert result.n_iter == 3


def test_result_from_pygimli_missing_chi2_and_iter_default():
    class Inv:
        pass

    cfg = InversionConfig(
        method="mt", dimension="1d", backend="pygimli",
        data={"freqs": [1.0], "rho_a": [100.0]},
    )
    em_data = EMData.coerce({"freqs": [1.0], "rho_a": [100.0]})
    recovered = StartingModel([100.0, 200.0], [50.0])
    result = pygimli_mod._result_from_pygimli(
        cfg, em_data, recovered, np.array([100.0]), 0.1,
        Inv(), fop="FOP", station_index=None, extra={"mode": "mt"},
    )
    assert np.isnan(result.objective)
    assert result.n_iter == 0


# ---------------------------------------------------------------------------
# _layer_centers / _station_data / _row / _station_names / _station_x
# ---------------------------------------------------------------------------


def test_layer_centers_normal():
    centers = pygimli_mod._layer_centers(np.array([10.0, 20.0]))
    assert centers.size == 3
    assert centers[0] == 5.0


def test_layer_centers_empty_thicknesses():
    centers = pygimli_mod._layer_centers(np.array([]))
    assert centers.size == 1
    assert centers[0] == 0.5


def test_station_data_builds_single_station_emdata():
    em_data = EMData(
        method="mt",
        frequencies=[1.0, 10.0],
        rho_a=[[80.0, 90.0], [70.0, 60.0]],
        phase=[[40.0, 41.0], [42.0, 43.0]],
        errors=[[1.0, 1.0], [2.0, 2.0]],
        station_names=["A", "B"],
        station_x=[0.0, 100.0],
    )
    single = pygimli_mod._station_data(em_data, 1)
    assert single.station_names == ["B"]
    np.testing.assert_allclose(single.rho_a, [70.0, 60.0])
    np.testing.assert_allclose(single.phase, [42.0, 43.0])
    np.testing.assert_allclose(single.station_x, [100.0])


def test_row_1d_and_2d_branches():
    arr_1d = np.array([1.0, 2.0, 3.0])
    np.testing.assert_allclose(pygimli_mod._row(arr_1d, 1), arr_1d)
    arr_2d = np.array([[1.0, 2.0], [3.0, 4.0]])
    np.testing.assert_allclose(pygimli_mod._row(arr_2d, 1), [3.0, 4.0])


def test_station_names_defaults_and_explicit():
    em_data_default = EMData(method="mt", frequencies=[1.0], rho_a=[[80.0], [90.0]])
    assert pygimli_mod._station_names(em_data_default, 2) == ["S000", "S001"]
    em_data_named = EMData(
        method="mt", frequencies=[1.0], rho_a=[[80.0], [90.0]],
        station_names=["X1", "X2"],
    )
    assert pygimli_mod._station_names(em_data_named, 2) == ["X1", "X2"]


def test_station_x_defaults_and_explicit():
    em_data_default = EMData(method="mt", frequencies=[1.0], rho_a=[[80.0], [90.0]])
    np.testing.assert_allclose(
        pygimli_mod._station_x(em_data_default, 2), [0.0, 1.0]
    )
    em_data_x = EMData(
        method="mt", frequencies=[1.0], rho_a=[[80.0], [90.0]], station_x=[5.0, 15.0]
    )
    np.testing.assert_allclose(pygimli_mod._station_x(em_data_x, 2), [5.0, 15.0])


# ---------------------------------------------------------------------------
# _load_pygimli -- real absence
# ---------------------------------------------------------------------------


def test_load_pygimli_raises_import_error_when_absent():
    with pytest.raises(ImportError, match="pyGIMLi"):
        pygimli_mod._load_pygimli()
