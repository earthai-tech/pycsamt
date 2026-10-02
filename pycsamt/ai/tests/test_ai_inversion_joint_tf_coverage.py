# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TensorFlow-backend coverage for pycsamt.ai.inversion.joint.JointInverter.

Requires TensorFlow to be installed and its native runtime to actually
load — run this file in an environment that has it (e.g. the
``base-attentiv``-adjacent ``pycsamt-tf-verify`` env used for the other
TF-backend coverage files in this batch); it is skipped everywhere
else. Same subprocess-probe rationale as
``pycsamt/ai/tests/test_pinn_ops_tf.py`` / ``test_ai_imputer_tf_coverage.py``.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest


def _tf_importable() -> bool:
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


def _synthetic_joint_data(n=60, n_features_list=(12, 8), n_layers=3, seed=0):
    rng = np.random.default_rng(seed)
    X_list = [
        rng.standard_normal((n, nf)).astype(np.float32)
        for nf in n_features_list
    ]
    n_out = 2 * n_layers - 1
    y = rng.uniform(1.0, 3.0, size=(n, n_out)).astype(np.float32)
    y[:, n_layers:] = rng.uniform(10.0, 200.0, size=(n, n_layers - 1))
    return X_list, y


def _new_inverter(**kwargs):
    from pycsamt.ai.inversion.joint import JointInverter

    return JointInverter(
        n_features_list=(12, 8),
        n_layers=3,
        hidden_dim=16,
        growth_rate=8,
        n_dense_layers=2,
        **kwargs,
    )


def test_fit_selects_tensorflow_backend():
    X_list, y = _synthetic_joint_data()
    inv = _new_inverter()
    inv.fit(X_list, y, epochs=2, batch_size=16, verbose=False)
    assert inv._backend_name == "tensorflow"
    assert inv._is_fitted is True
    assert "train_loss" in inv._history
    assert "val_loss" in inv._history


def test_predict_shape_tensorflow():
    X_list, y = _synthetic_joint_data()
    inv = _new_inverter()
    inv.fit(X_list, y, epochs=2, batch_size=16, verbose=False)
    y_pred = inv.predict(X_list)
    assert y_pred.shape == (60, 5)
    assert np.all(np.isfinite(y_pred))


def test_save_load_round_trip_tensorflow(tmp_path):
    from pycsamt.ai.inversion.joint import JointInverter

    X_list, y = _synthetic_joint_data()
    inv = _new_inverter()
    inv.fit(X_list, y, epochs=2, batch_size=16, verbose=False)
    y_before = inv.predict(X_list)

    path = tmp_path / "joint_tf.npz"
    inv.save(path)
    loaded = JointInverter.load(path)

    assert loaded._backend_name == "tensorflow"
    y_after = loaded.predict(X_list)
    np.testing.assert_allclose(y_before, y_after, rtol=1e-4, atol=1e-5)
