"""Real-fit tests for EMInverter2D's staged spatial-regularization loss.

Unlike test_ai_inversion_api_contracts.py, these tests actually train
small networks to prove the new lambda_x/lambda_z/lambda_tv wiring
works end to end, not just that it imports.
"""

from __future__ import annotations

import numpy as np
import pytest
from unittest.mock import patch

torch = pytest.importorskip("torch")

from pycsamt.ai.inversion.inv2d import EMInverter2D, _staged_torch_loss  # noqa: E402


def _toy_data(seed=0, n=12, n_freq=8, n_sta=6, n_depth=8):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, 2, n_freq, n_sta)).astype(np.float32)
    y = rng.normal(size=(n, n_depth, n_sta)).astype(np.float32)
    return X, y


def test_staged_loss_matches_plain_mse_at_zero_weights():
    pred = torch.randn(2, 1, 6, 5, requires_grad=True)
    target = torch.randn(2, 1, 6, 5)
    plain = torch.nn.functional.mse_loss(pred, target)
    staged = _staged_torch_loss(pred, target)
    assert torch.allclose(plain, staged)


def test_staged_loss_changes_with_nonzero_weights_and_is_differentiable():
    pred = torch.randn(2, 1, 6, 5, requires_grad=True)
    target = torch.randn(2, 1, 6, 5)
    baseline = _staged_torch_loss(pred, target)
    staged = _staged_torch_loss(
        pred, target, lambda_x=0.1, lambda_z=0.1, lambda_tv=0.05
    )
    assert not torch.allclose(baseline, staged)
    staged.backward()
    assert torch.isfinite(pred.grad).all()


def test_fit_default_lambdas_record_zero_weights():
    X, y = _toy_data()
    inv = EMInverter2D(n_components=2, n_depth=8, n_stations=6, n_freqs=8)
    inv.fit(X, y, epochs=3, batch_size=4, patience=10, verbose=False, seed=0)
    assert inv._meta["loss_weights"] == {
        "lambda_x": 0.0,
        "lambda_z": 0.0,
        "lambda_tv": 0.0,
    }
    assert np.isfinite(inv._meta["best_val_loss"])


def test_fit_with_staged_weights_trains_and_predicts():
    X, y = _toy_data()
    inv = EMInverter2D(n_components=2, n_depth=8, n_stations=6, n_freqs=8)
    inv.fit(
        X,
        y,
        epochs=3,
        batch_size=4,
        patience=10,
        verbose=False,
        seed=0,
        lambda_x=0.1,
        lambda_z=0.1,
        lambda_tv=0.05,
    )
    assert inv._is_fitted
    assert inv._meta["loss_weights"] == {
        "lambda_x": 0.1,
        "lambda_z": 0.1,
        "lambda_tv": 0.05,
    }
    pred = inv.predict(X[:2])
    assert pred.shape == (2, 8, 6)
    assert np.all(np.isfinite(pred))


def test_fit_rejects_staged_weights_on_tensorflow_backend():
    X, y = _toy_data()
    inv = EMInverter2D(n_components=2, n_depth=8, n_stations=6, n_freqs=8)
    with patch(
        "pycsamt.ai.inversion.inv2d.active_backend", return_value="tensorflow"
    ):
        with pytest.raises(NotImplementedError, match="PyTorch backend"):
            inv.fit(X, y, epochs=1, lambda_x=0.1, verbose=False)
