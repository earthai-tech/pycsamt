# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.backends._torch.TorchBackend.

Architecture forward-pass correctness is already covered by
``pycsamt/ai/tests/test_nets.py`` — these tests exercise the backend
*wrapper* itself: dispatch, device resolution, the train/predict loop,
weight round-tripping, and loss-function selection.
"""

from __future__ import annotations

import importlib.util
from unittest import mock

import numpy as np
import pytest

_HAS_TORCH = importlib.util.find_spec("torch") is not None

pytestmark = pytest.mark.skipif(not _HAS_TORCH, reason="PyTorch not installed")


@pytest.fixture
def backend():
    from pycsamt.backends._torch import TorchBackend

    return TorchBackend()


# ─── availability / require ─────────────────────────────────────────────────


def test_is_available_true(backend):
    assert backend.is_available() is True


def test_require_does_not_raise(backend):
    backend.require()  # torch installed → no exception


def test_require_raises_when_unavailable():
    from pycsamt.backends._torch import TorchBackend

    with mock.patch.object(TorchBackend, "is_available", return_value=False):
        with pytest.raises(ImportError, match="pip install torch"):
            TorchBackend.require()


# ─── resolve_device ──────────────────────────────────────────────────────────


def test_resolve_device_explicit_passthrough(backend):
    assert backend.resolve_device("cuda:1") == "cuda:1"


def test_resolve_device_cpu_when_nothing_available(backend):
    import torch

    with (
        mock.patch.object(torch.cuda, "is_available", return_value=False),
        mock.patch.object(torch.backends.mps, "is_available", return_value=False),
    ):
        assert backend.resolve_device(None) == "cpu"


def test_resolve_device_cuda_preferred(backend):
    import torch

    with mock.patch.object(torch.cuda, "is_available", return_value=True):
        assert backend.resolve_device(None) == "cuda"


def test_resolve_device_mps_when_no_cuda(backend):
    import torch

    if not hasattr(torch.backends, "mps"):
        pytest.skip("torch build has no mps backend attribute")
    with (
        mock.patch.object(torch.cuda, "is_available", return_value=False),
        mock.patch.object(torch.backends.mps, "is_available", return_value=True),
    ):
        assert backend.resolve_device(None) == "mps"


def test_resolve_device_import_error_falls_back_to_cpu(backend):
    """If `import torch` itself fails inside resolve_device, it must be
    swallowed and 'cpu' returned rather than propagating."""
    import sys

    with mock.patch.dict(sys.modules, {"torch": None}):
        assert backend.resolve_device(None) == "cpu"


# ─── build dispatch ──────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "arch,spec_extra",
    [
        ("cnn1d", {"n_features": 20, "n_out": 5}),
        ("cnn", {"n_features": 20, "n_out": 5}),
        ("resnet1d", {"n_features": 20, "n_out": 5}),
        ("resnet", {"n_features": 20, "n_out": 5}),
        ("fcn1d", {"n_features": 20, "n_out": 5}),
        ("fcn", {"n_features": 20, "n_out": 5}),
        ("unet2d", {"n_in": 2, "n_out": 1}),
        ("unet", {"n_in": 2, "n_out": 1}),
        ("drcnn", {"n_features_list": (10, 8), "n_out": 3}),
        ("CNN1D", {"n_features": 20, "n_out": 5}),  # case-insensitive
    ],
)
def test_build_dispatches_to_correct_architecture(backend, arch, spec_extra):
    import torch.nn as nn

    spec = {"arch": arch, **spec_extra}
    model = backend.build(spec)
    assert isinstance(model, nn.Module)


def test_build_unknown_arch_raises(backend):
    with pytest.raises(ValueError, match="Unknown arch"):
        backend.build({"arch": "not_a_real_arch"})


# ─── train / predict / weights round-trip ──────────────────────────────────


@pytest.fixture
def small_data():
    rng = np.random.default_rng(0)
    X_train = rng.standard_normal((16, 6)).astype(np.float32)
    y_train = rng.standard_normal((16, 2)).astype(np.float32)
    X_val = rng.standard_normal((8, 6)).astype(np.float32)
    y_val = rng.standard_normal((8, 2)).astype(np.float32)
    return X_train, y_train, X_val, y_val


def test_train_runs_and_returns_history(backend, small_data):
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "cnn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    hist = backend.train(
        model,
        X_train,
        y_train,
        X_val,
        y_val,
        epochs=2,
        batch_size=4,
        patience=1,
        verbose=False,
        device="cpu",
    )
    assert set(hist) == {"train_loss", "val_loss"}
    assert len(hist["train_loss"]) == len(hist["val_loss"])
    assert len(hist["train_loss"]) >= 1


def test_train_early_stopping_triggers(backend, small_data):
    """patience=0 must stop after the first non-improving epoch."""
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "fcn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    hist = backend.train(
        model,
        X_train,
        y_train,
        X_val,
        y_val,
        epochs=50,
        batch_size=4,
        patience=0,
        verbose=False,
        device="cpu",
    )
    assert len(hist["train_loss"]) < 50


def test_train_with_grad_clip_and_verbose(backend, small_data, capsys):
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "cnn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    hist = backend.train(
        model,
        X_train,
        y_train,
        X_val,
        y_val,
        epochs=1,
        batch_size=4,
        patience=5,
        grad_clip=1.0,
        verbose=True,
        device="cpu",
    )
    assert len(hist["train_loss"]) == 1
    # verbose=True now renders an animated tqdm bar (see
    # pycsamt.utils.progress.ProgressBar), which — like tqdm itself —
    # writes to stderr, not stdout, so stdout stays clean for piping.
    captured = capsys.readouterr()
    assert "torch" in captured.err


def test_train_zero_epochs_leaves_best_state_none(backend, small_data):
    """epochs=0 never enters the training loop: best_state stays None and
    the `if best_state is not None` guard must be skipped without error."""
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "cnn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    hist = backend.train(
        model,
        X_train,
        y_train,
        X_val,
        y_val,
        epochs=0,
        batch_size=4,
        patience=5,
        verbose=False,
        device="cpu",
    )
    assert hist == {"train_loss": [], "val_loss": []}


def test_train_early_stop_prints_message(backend, small_data, capsys):
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "cnn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    backend.train(
        model,
        X_train,
        y_train,
        X_val,
        y_val,
        epochs=10,
        batch_size=4,
        patience=0,
        verbose=True,
        device="cpu",
    )
    out = capsys.readouterr().out
    assert "Early stop at epoch" in out


def test_train_over_many_epochs_covers_non_improving_step(backend, small_data):
    """With enough epochs on noisy random data, at least one epoch fails to
    beat the running best val loss — exercising the no_improve += 1 path."""
    import torch

    torch.manual_seed(0)
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "fcn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    hist = backend.train(
        model,
        X_train,
        y_train,
        X_val,
        y_val,
        epochs=25,
        batch_size=4,
        patience=100,
        verbose=False,
        device="cpu",
    )
    assert len(hist["val_loss"]) == 25


def test_predict_shape_and_batching(backend, small_data):
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "cnn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    preds = backend.predict(model, X_val, batch_size=3)
    assert preds.shape == (8, 2)
    assert preds.dtype == np.float32


def test_get_set_weights_roundtrip(backend):
    model_a = backend.build(
        {"arch": "cnn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    model_b = backend.build(
        {"arch": "cnn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    weights = backend.get_weights(model_a)
    assert all(isinstance(v, np.ndarray) for v in weights.values())

    backend.set_weights(model_b, weights)
    weights_b = backend.get_weights(model_b)

    assert set(weights) == set(weights_b)
    for k in weights:
        np.testing.assert_array_equal(weights[k], weights_b[k])


def test_predict_matches_after_weight_transfer(backend):
    rng = np.random.default_rng(1)
    X = rng.standard_normal((5, 6)).astype(np.float32)

    model_a = backend.build(
        {"arch": "fcn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    model_b = backend.build(
        {"arch": "fcn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    backend.set_weights(model_b, backend.get_weights(model_a))

    pred_a = backend.predict(model_a, X)
    pred_b = backend.predict(model_b, X)
    np.testing.assert_allclose(pred_a, pred_b, rtol=1e-5, atol=1e-6)


# ─── loss-function selection ────────────────────────────────────────────────


def test_get_loss_fn_mse(backend):
    import torch.nn as nn

    fn = backend._get_loss_fn("mse")
    assert isinstance(fn, nn.MSELoss)


def test_get_loss_fn_masked_mse(backend):
    import torch

    fn = backend._get_loss_fn("masked_mse")
    pred = torch.tensor([[1.0, 2.0]])
    target = torch.tensor([[1.0, float("nan")]])
    loss = fn(pred, target)
    assert torch.isfinite(loss)


def test_get_loss_fn_cross_entropy(backend):
    import torch.nn as nn

    fn = backend._get_loss_fn("cross_entropy")
    assert isinstance(fn, nn.CrossEntropyLoss)


def test_get_loss_fn_unknown_raises(backend):
    with pytest.raises(ValueError, match="Unknown loss"):
        backend._get_loss_fn("not_a_loss")


def test_train_masked_mse_loss(backend, small_data):
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "fcn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    hist = backend.train(
        model,
        X_train,
        y_train,
        X_val,
        y_val,
        epochs=1,
        batch_size=4,
        patience=1,
        loss="masked_mse",
        verbose=False,
        device="cpu",
    )
    assert len(hist["train_loss"]) == 1


# ─── repr (inherited from NeuralBackend) ────────────────────────────────────


def test_repr_available(backend):
    assert "available" in repr(backend)
