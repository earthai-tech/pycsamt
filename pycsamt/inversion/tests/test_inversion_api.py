# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from pycsamt.forward import MT1DForward, TEM1DForward, LayeredModel
from pycsamt.inversion import (
    EMData,
    InversionConfig,
    InversionWorkflow,
    StartingModel,
    available_backends,
    export,
    plot,
)
from pycsamt.inversion.mt import MT1DInversion
from pycsamt.inversion.mt import MT2DInversion, MT3DInversion
from pycsamt.inversion.amt import AMT1DInversion, AMT2DInversion
from pycsamt.inversion.emap import EMAP2DInversion
from pycsamt.inversion.tdem import TDEM1DInversion, TDEM2DInversion
from pycsamt.inversion.backends import simpeg as simpeg_backend
from pycsamt.inversion.backends import pygimli as pygimli_backend


def test_config_template_roundtrip(tmp_path):
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="builtin",
        n_layers=3,
        workdir=str(tmp_path / "run"),
    )
    path = cfg.write_template(tmp_path / "inv_config.py")
    loaded = InversionConfig.from_file(path)
    assert loaded.method == "mt"
    assert loaded.dimension == "1d"
    assert loaded.n_layers == 3


def test_emdata_coerce_aliases():
    data = EMData.coerce(
        {
            "method": "AMT",
            "freqs": [1.0, 10.0],
            "rho_a": [100.0, 120.0],
            "phase": [45.0, 47.0],
            "stations": ["S1"],
        }
    )
    assert data.method == "amt"
    assert data.has_mt_response
    assert data.n_samples == 2
    assert data.station_names == ["S1"]


def test_backend_registry():
    names = available_backends()
    assert "builtin" in names
    assert "occam2d" in names
    assert "simpeg" in names


def test_simpeg_backend_missing_dependency_is_clear():
    if importlib.util.find_spec("simpeg") is not None:
        pytest.skip("SimPEG is installed; missing-dependency path not applicable.")
    cfg = InversionConfig(
        method="mt",
        dimension="1d",
        backend="simpeg",
        data={"freqs": [1.0], "rho_a": [100.0], "phase": [45.0]},
    )
    with pytest.raises(ImportError, match="SimPEG"):
        InversionWorkflow(cfg).run()


def test_simpeg_observation_packing_order_and_errors():
    data = EMData.coerce({
        "freqs": [1.0, 10.0],
        "rho_a": [100.0, 200.0],
        "phase": [45.0, 50.0],
    })
    cfg = InversionConfig(error_floor=0.05, phase_error=2.0)
    values, errors = simpeg_backend._pack_nsem_observations(data, cfg)
    np.testing.assert_allclose(values, [100.0, 45.0, 200.0, 50.0])
    np.testing.assert_allclose(errors, [5.0, 2.0, 10.0, 2.0])


def test_simpeg_3d_observation_packing_and_station_locations():
    data = EMData.coerce({
        "freqs": [1.0, 10.0],
        "rho_a": [[100.0, 200.0], [110.0, 210.0]],
        "phase": [[45.0, 50.0], [46.0, 51.0]],
        "station_x": [0.0, 500.0],
        "metadata": {"station_y": [10.0, 20.0]},
    })
    cfg = InversionConfig(error_floor=0.05, phase_error=2.0)
    values, errors = simpeg_backend._pack_nsem_observations(data, cfg)
    np.testing.assert_allclose(
        values,
        [100.0, 110.0, 45.0, 46.0, 200.0, 210.0, 50.0, 51.0],
    )
    np.testing.assert_allclose(errors, [5.0, 5.5, 2.0, 2.0, 10.0, 10.5, 2.0, 2.0])
    locations = simpeg_backend._station_locations(data)
    np.testing.assert_allclose(locations, [[0.0, 10.0, 0.0], [500.0, 20.0, 0.0]])


def test_simpeg_3d_model_reshape_helper():
    class Mesh:
        shape_cells = (2, 2, 2)

    log_sigma = np.log(np.arange(1.0, 9.0))
    rho = simpeg_backend._rho_3d_from_log_sigma(log_sigma, Mesh())
    assert rho.shape == (2, 2, 2)
    np.testing.assert_allclose(
        np.sort(rho.reshape(-1)),
        np.sort(1.0 / np.arange(1.0, 9.0)),
    )


def test_pygimli_backend_missing_dependency_is_clear():
    if importlib.util.find_spec("pygimli") is not None:
        pytest.skip("pyGIMLi is installed; missing-dependency path not applicable.")
    cfg = InversionConfig(
        method="tdem",
        dimension="1d",
        backend="pygimli",
        data={"times": [1e-4], "values": [1e-9]},
    )
    with pytest.raises(ImportError, match="pyGIMLi"):
        InversionWorkflow(cfg).run()


def test_pygimli_observation_packing_helpers():
    data = EMData.coerce({
        "freqs": [1.0, 10.0],
        "rho_a": [100.0, 200.0],
        "phase": [45.0, 50.0],
    })
    values = pygimli_backend._pack_mt_observations(data)
    errors = pygimli_backend._pack_mt_errors(
        data,
        InversionConfig(error_floor=0.05, phase_error=2.0),
    )
    np.testing.assert_allclose(values, [100.0, 200.0, 45.0, 50.0])
    np.testing.assert_allclose(errors[:2], [0.05, 0.05])
    assert np.all(errors[2:] > 0)


def test_expected_phase_tree_modules_exist():
    import importlib

    modules = [
        "pycsamt.inversion.mt.inv1d",
        "pycsamt.inversion.mt.inv2d",
        "pycsamt.inversion.mt.inv3d",
        "pycsamt.inversion.amt.inv1d",
        "pycsamt.inversion.amt.inv2d",
        "pycsamt.inversion.emap.inv2d",
        "pycsamt.inversion.tdem.inv1d",
        "pycsamt.inversion.tdem.inv2d",
    ]
    for module in modules:
        assert importlib.import_module(module)


def test_method_wrapper_defaults_do_not_duplicate_constructor_logic():
    wrappers = [
        (MT1DInversion, "mt", "1d", "builtin"),
        (MT2DInversion, "mt", "2d", "builtin"),
        (MT3DInversion, "mt", "3d", "modem"),
        (AMT1DInversion, "amt", "1d", "builtin"),
        (AMT2DInversion, "amt", "2d", "builtin"),
        (EMAP2DInversion, "emap", "2d", "builtin"),
        (TDEM1DInversion, "tdem", "1d", "builtin"),
        (TDEM2DInversion, "tdem", "2d", "builtin"),
    ]
    for cls, method, dimension, backend in wrappers:
        inv = cls(max_iter=1)
        assert inv.config.method == method
        assert inv.config.dimension == dimension
        assert inv.config.backend == backend


def test_builtin_1d_inversion_smoke(tmp_path):
    freqs = np.logspace(-2, 2, 12)
    true_model = LayeredModel([80.0, 20.0, 500.0], [250.0, 900.0])
    response = MT1DForward(freqs=freqs).run(true_model)
    data = {
        "freqs": freqs,
        "rho_a": response.rho_a,
        "phase": response.phase,
    }
    start = StartingModel([100.0, 50.0, 300.0], [200.0, 1000.0])
    cfg = InversionConfig(
        data=data,
        n_layers=3,
        starting_model=start,
        max_iter=25,
        workdir=str(tmp_path / "builtin"),
    )
    result = InversionWorkflow(cfg).run()
    assert result.backend == "builtin"
    assert result.model.n_layers == 3
    assert np.isfinite(result.rms)
    assert result.to_resistivity_model().rho_2d.shape == (3, 1)


def test_mt_convenience_wrapper(tmp_path):
    freqs = np.logspace(-1, 1, 5)
    response = MT1DForward(freqs=freqs).run(
        LayeredModel([100.0, 300.0], [500.0])
    )
    inv = MT1DInversion(
        {"freqs": freqs, "rho_a": response.rho_a, "phase": response.phase},
        n_layers=2,
        max_iter=5,
        workdir=str(tmp_path / "mt"),
    )
    result = inv.run()
    assert result.method == "mt"
    assert result.dimension == "1d"


def test_builtin_tdem_1d_inversion_smoke(tmp_path):
    times = np.logspace(-5, -3, 4)
    true_model = LayeredModel([80.0, 300.0], [120.0])
    fwd_options = {"loop_radius": 25.0, "n_freqs": 8, "n_lam": 12}
    response = TEM1DForward(times=times, **fwd_options).run(true_model)
    inv = TDEM1DInversion(
        {
            "times": times,
            "values": response.dBz_dt,
        },
        n_layers=2,
        starting_model=StartingModel([100.0, 250.0], [150.0]),
        max_iter=2,
        backend_options=fwd_options,
        workdir=str(tmp_path / "tdem"),
    )
    result = inv.run()
    assert result.method == "tdem"
    assert result.dimension == "1d"
    assert result.backend == "builtin"
    assert result.model.n_layers == 2
    assert np.isfinite(result.rms)


def test_builtin_2d_stitched_profile_and_exports(tmp_path):
    freqs = np.logspace(-2, 2, 10)
    models = [
        LayeredModel([80.0, 20.0, 500.0], [250.0, 900.0]),
        LayeredModel([100.0, 35.0, 600.0], [300.0, 850.0]),
        LayeredModel([120.0, 60.0, 700.0], [350.0, 800.0]),
    ]
    responses = [MT1DForward(freqs=freqs).run(model) for model in models]
    data = {
        "method": "mt",
        "freqs": freqs,
        "rho_a": np.vstack([resp.rho_a for resp in responses]),
        "phase": np.vstack([resp.phase for resp in responses]),
        "station_names": ["S1", "S2", "S3"],
        "station_x": [0.0, 500.0, 1000.0],
    }
    cfg = InversionConfig(
        method="mt",
        dimension="2d",
        backend="builtin",
        data=data,
        n_layers=3,
        starting_model=StartingModel([100.0, 50.0, 500.0], [300.0, 900.0]),
        max_iter=15,
        workdir=str(tmp_path / "profile"),
    )
    result = InversionWorkflow(cfg).run()
    model = result.to_resistivity_model()
    assert result.dimension == "2d"
    assert model.rho_2d.shape == (3, 3)
    assert model.station_names == ["S1", "S2", "S3"]
    assert np.isfinite(result.rms)

    csv_path = export.to_csv(result, tmp_path / "profile.csv")
    npz_path = export.to_npz(result, tmp_path / "profile.npz")
    assert csv_path.exists()
    assert npz_path.exists()

    ax_model = plot.plot_model(result, section="compact", colorbar=False)
    ax_rms = plot.plot_rms(result)
    assert ax_model.get_ylabel() == "Depth (m)"
    assert ax_rms.get_ylabel() == "Weighted RMS"
