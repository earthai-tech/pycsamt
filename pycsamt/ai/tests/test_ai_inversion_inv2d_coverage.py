# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Coverage tests for pycsamt.ai.inversion.inv2d.EMInverter2D.

test_ai_inversion_api_contracts.py only exercises the no-backend
interface (construction, ``_get_params``, predict-before-fit).  This
file drives a full PyTorch ``fit``/``predict`` cycle on tiny synthetic
2-D panels: the staged spatial-regularization loss
(``lambda_x``/``lambda_z``/``lambda_tv``), the TensorFlow
staged-regularization guard (monkeypatched, no real TensorFlow import),
``_resolve_channels``, and the save/load round trip.
"""

from __future__ import annotations

import numpy as np
import pytest


def _has_backend():
    from pycsamt.backends import get_backend

    return get_backend() != "none"


pytestmark = pytest.mark.skipif(not _has_backend(), reason="no DL backend available")


def _synthetic_2d_data(
    n=10, n_components=2, n_depth=8, n_stations=8, n_freqs=8, seed=0
):
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, n_components, n_freqs, n_stations)).astype(
        np.float32
    )
    y = rng.uniform(1.0, 3.0, size=(n, n_depth, n_stations)).astype(
        np.float32
    )
    return X, y


def _new_inv(**kwargs):
    from pycsamt.ai.inversion.inv2d import EMInverter2D

    defaults = dict(
        n_components=2,
        n_depth=8,
        n_stations=8,
        n_freqs=8,
        channels=(4, 8),
    )
    defaults.update(kwargs)
    return EMInverter2D(**defaults)


def test_resolve_channels_explicit_vs_auto():
    inv_explicit = _new_inv(channels=(4, 8, 16))
    assert inv_explicit._channels == (4, 8, 16)

    inv_auto = _new_inv(channels=None, n_freqs=32, n_stations=32)
    # auto channel spec has n_stages encoder widths + 1 bridge width
    assert len(inv_auto._channels) >= 2

    inv_auto_depth = _new_inv(channels=None, unet_depth=2, n_freqs=32, n_stations=32)
    assert len(inv_auto_depth._channels) == 3  # 2 stages + 1 bridge


def test_fit_marks_fitted_and_history():
    X, y = _synthetic_2d_data()
    inv = _new_inv()
    inv.fit(X, y, epochs=2, batch_size=4, verbose=False)
    assert inv._is_fitted
    assert "train_loss" in inv._history
    assert "val_loss" in inv._history
    assert "best_val_loss" in inv._meta
    assert inv._backend_name == "torch"
    assert "fitted" in repr(inv)


def test_predict_shape_and_log_linear_consistency():
    X, y = _synthetic_2d_data()
    inv = _new_inv()
    inv.fit(X, y, epochs=2, batch_size=4, verbose=False)

    y_log = inv.predict(X, as_log_rho=True)
    y_lin = inv.predict(X, as_log_rho=False)
    assert y_log.shape == (10, 8, 8)
    np.testing.assert_allclose(10.0**y_log, y_lin, rtol=1e-4)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"lambda_x": 0.1},
        {"lambda_z": 0.1},
        {"lambda_tv": 0.1},
        {"lambda_x": 0.05, "lambda_z": 0.05, "lambda_tv": 0.05},
    ],
)
def test_fit_with_staged_regularization_terms(kwargs):
    X, y = _synthetic_2d_data(n=8)
    inv = _new_inv()
    inv.fit(X, y, epochs=2, batch_size=4, verbose=False, **kwargs)
    assert inv._is_fitted
    assert inv._meta["loss_weights"] == {
        "lambda_x": kwargs.get("lambda_x", 0.0),
        "lambda_z": kwargs.get("lambda_z", 0.0),
        "lambda_tv": kwargs.get("lambda_tv", 0.0),
    }


def test_fit_grad_clip_none():
    X, y = _synthetic_2d_data(n=8)
    inv = _new_inv()
    inv.fit(X, y, epochs=2, batch_size=4, grad_clip=None, verbose=False)
    assert inv._is_fitted


def test_tensorflow_backend_rejects_staged_regularization(monkeypatch):
    """Guards the TF path without importing real TensorFlow (in-process
    TF import can crash this environment -- see fix_tensorflow_broken_dll
    memory note). Only the early guard-clause branch is exercised."""
    import pycsamt.ai.inversion.inv2d as inv2d_mod

    monkeypatch.setattr(inv2d_mod, "active_backend", lambda: "tensorflow")

    X, y = _synthetic_2d_data(n=4)
    inv = _new_inv()
    with pytest.raises(NotImplementedError, match="PyTorch backend"):
        inv.fit(X, y, epochs=1, lambda_x=0.1, verbose=False)


def test_save_load_round_trip(tmp_path):
    from pycsamt.ai.inversion.inv2d import EMInverter2D

    X, y = _synthetic_2d_data(n=8)
    inv = _new_inv()
    inv.fit(X, y, epochs=2, batch_size=4, verbose=False)
    y_before = inv.predict(X)

    path = tmp_path / "inv2d.npz"
    inv.save(path)
    loaded = EMInverter2D.load(path)

    y_after = loaded.predict(X)
    np.testing.assert_allclose(y_before, y_after, rtol=1e-4, atol=1e-5)
    assert loaded._channels == inv._channels


def test_predict_before_fit_raises():
    inv = _new_inv()
    with pytest.raises(RuntimeError, match=r"fit\(\) before predict"):
        inv.predict(np.zeros((1, 2, 8, 8), dtype=np.float32))
