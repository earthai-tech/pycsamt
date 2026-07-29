# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.backends._tensorflow.TensorFlowBackend.

Requires TensorFlow to be installed — run this file in an environment
that has it (e.g. the ``base-attentiv`` conda env), it is skipped
everywhere else.
"""

from __future__ import annotations

import subprocess
import sys
from unittest import mock

import numpy as np
import pytest


def _tf_importable() -> bool:
    """Probe TensorFlow importability out-of-process.

    `find_spec` alone is not enough: the package can be installed (spec
    found) while the native runtime fails to load (e.g. a broken DLL) —
    that must skip, not fail, this file. Importing it in-process to check
    is *not* safe either: on a broken install, `import tensorflow` can
    trigger a native access violation that crashes the whole pytest
    process outright (observed with a DLL conflict against an
    already-loaded torch/numpy here) rather than raising a catchable
    Python exception. Running the probe in a subprocess isolates any such
    crash from the test runner.
    """
    try:
        result = subprocess.run(
            [sys.executable, "-c", "import tensorflow"],
            capture_output=True,
            timeout=60,
        )
        return result.returncode == 0
    except Exception:
        return False


_HAS_TF = _tf_importable()

pytestmark = pytest.mark.skipif(
    not _HAS_TF, reason="TensorFlow not installed or not importable"
)


@pytest.fixture
def backend():
    from pycsamt.backends._tensorflow import TensorFlowBackend

    return TensorFlowBackend()


# ─── availability / require ─────────────────────────────────────────────────


def test_is_available_true(backend):
    assert backend.is_available() is True


def test_require_does_not_raise(backend):
    backend.require()


def test_require_raises_when_unavailable():
    from pycsamt.backends._tensorflow import TensorFlowBackend

    with mock.patch.object(TensorFlowBackend, "is_available", return_value=False):
        with pytest.raises(ImportError, match="pip install tensorflow"):
            TensorFlowBackend.require()


# ─── resolve_device ──────────────────────────────────────────────────────────


def test_resolve_device_explicit_passthrough(backend):
    assert backend.resolve_device("/GPU:1") == "/GPU:1"


def test_resolve_device_cpu_when_no_gpu(backend):
    import tensorflow as tf

    with mock.patch.object(tf.config, "list_physical_devices", return_value=[]):
        assert backend.resolve_device(None) == "/CPU:0"


def test_resolve_device_gpu_when_available(backend):
    import tensorflow as tf

    with mock.patch.object(tf.config, "list_physical_devices", return_value=["GPU:0"]):
        assert backend.resolve_device(None) == "/GPU:0"


def test_resolve_device_import_error_falls_back_to_cpu(backend):
    import sys

    with mock.patch.dict(sys.modules, {"tensorflow": None}):
        assert backend.resolve_device(None) == "/CPU:0"


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
        (
            "unet2d",
            {"n_in": 2, "n_out": 1, "input_shape": (16, 16, 2)},
        ),
        ("unet", {"n_in": 2, "n_out": 1, "input_shape": (16, 16, 2)}),
        ("drcnn", {"n_features_list": (10, 8), "n_out": 3}),
        ("CNN1D", {"n_features": 20, "n_out": 5}),  # case-insensitive
    ],
)
def test_build_dispatches_to_correct_architecture(backend, arch, spec_extra):
    import tensorflow as tf

    spec = {"arch": arch, **spec_extra}
    model = backend.build(spec)
    assert isinstance(model, tf.keras.Model)


def test_build_unknown_arch_raises(backend):
    with pytest.raises(ValueError, match="Unknown arch"):
        backend.build({"arch": "not_a_real_arch"})


def test_build_resnet1d_multi_stage_hits_projection_and_identity_shortcut(
    backend,
):
    """Default channels=(64,128,256) exercise both the identity shortcut
    (first block, matching channel count) and the projected shortcut
    (stage transitions with stride=2 / channel change)."""
    model = backend.build({"arch": "resnet1d", "n_features": 20, "n_out": 4})
    out = model(np.zeros((2, 20), dtype=np.float32))
    assert out.shape == (2, 4)


def test_build_unet2d_default_input_shape(backend):
    """No explicit input_shape falls back to (None, None, n_in)."""
    model = backend.build({"arch": "unet2d", "n_in": 3, "n_out": 1})
    assert model.input_shape == (None, None, None, 3)


def test_build_drcnn_single_modality(backend):
    """len(encoded) == 1 skips the multi-modality Concatenate branch."""
    model = backend.build({"arch": "drcnn", "n_features_list": (10,), "n_out": 3})
    out = model(np.zeros((2, 10), dtype=np.float32))
    assert out.shape == (2, 3)


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
        device="/CPU:0",
    )
    assert set(hist) == {"train_loss", "val_loss"}
    assert len(hist["train_loss"]) == len(hist["val_loss"]) == 2


def test_train_verbose_true(backend, small_data):
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
        verbose=True,
        device="/CPU:0",
    )
    assert len(hist["train_loss"]) == 1


def test_train_default_device(backend, small_data):
    """device=None exercises resolve_device() inside train()."""
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "cnn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    hist = backend.train(model, X_train, y_train, X_val, y_val, epochs=1, batch_size=4)
    assert len(hist["train_loss"]) == 1


def test_predict_shape_and_batching(backend, small_data):
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "cnn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    preds = backend.predict(model, X_val, batch_size=3, device="/CPU:0")
    assert preds.shape == (8, 2)


def test_predict_default_device(backend, small_data):
    X_train, y_train, X_val, y_val = small_data
    model = backend.build(
        {"arch": "cnn1d", "n_features": 6, "n_out": 2, "channels": (4,)}
    )
    preds = backend.predict(model, X_val)
    assert preds.shape == (8, 2)


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
    assert backend._get_loss_fn("mse") == "mse"


def test_get_loss_fn_masked_mse_maps_to_mse(backend):
    assert backend._get_loss_fn("masked_mse") == "mse"


def test_get_loss_fn_cross_entropy(backend):
    assert backend._get_loss_fn("cross_entropy") == "sparse_categorical_crossentropy"


def test_get_loss_fn_passthrough_for_custom_loss(backend):
    assert backend._get_loss_fn("huber") == "huber"


# ─── repr (inherited from NeuralBackend) ────────────────────────────────────


def test_repr_available(backend):
    assert "available" in repr(backend)
