# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TensorFlow-backend coverage for pycsamt.ai.processing.imputer.EMImputer.

Requires TensorFlow to be installed and its native runtime to actually
load — run this file in an environment that has it (e.g. the
``base-attentiv``-adjacent ``pycsamt-tf-verify`` env used for
``_pinn_ops_tf.py``); it is skipped everywhere else.

Same subprocess-probe rationale as
``pycsamt/backends/tests/test_tensorflow_backend.py`` and
``pycsamt/ai/tests/test_pinn_ops_tf.py``: a broken native TF DLL can
crash the whole pytest process if probed with an in-process import.
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


def _imputer_X(n=20, n_comp=4, n_freqs=12, gap_frac=0.15, seed=4):
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, n_comp, n_freqs)).astype(np.float32)
    gaps = rng.random(X.shape) < gap_frac
    X[gaps] = np.nan
    return X, gaps


def test_fit_selects_tensorflow_backend():
    from pycsamt.ai.processing.imputer import EMImputer

    X, _ = _imputer_X()
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=2, batch_size=8, verbose=False)
    assert imp._backend_name == "tensorflow"
    assert imp._use_numpy is False
    assert set(imp.history_) == {"train_loss", "val_loss"}
    assert len(imp.history_["train_loss"]) == 2


def test_fit_tensorflow_verbose_prints(capsys):
    from pycsamt.ai.processing.imputer import EMImputer

    X, _ = _imputer_X(n=16, n_freqs=10)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=4, batch_size=8, verbose=True)
    out = capsys.readouterr().out
    assert "Epoch" in out


def test_fit_tensorflow_lr_plateau_reduction():
    from pycsamt.ai.processing.imputer import EMImputer

    X, _ = _imputer_X(n=24, n_freqs=8, gap_frac=0.3, seed=9)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=20, batch_size=8, lr=1e-2, verbose=False)
    assert len(imp.history_["val_loss"]) == 20


def test_transform_preserves_observed_and_fills_missing_tf():
    from pycsamt.ai.processing.imputer import EMImputer

    X, gaps = _imputer_X(n=20, n_freqs=12, gap_frac=0.2)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=3, batch_size=8, verbose=False)

    X_filled = imp.transform(X)
    assert X_filled.shape == X.shape
    assert np.all(np.isfinite(X_filled))
    np.testing.assert_array_equal(X_filled[~gaps], X[~gaps])
    assert not np.any(np.isnan(X_filled[gaps]))


def test_save_load_round_trip_tensorflow(tmp_path):
    from pycsamt.ai.processing.imputer import EMImputer

    X, _ = _imputer_X(n=16, n_freqs=10)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=2, batch_size=8, verbose=False)
    out_before = imp.transform(X)

    path = tmp_path / "imputer_tf.npz"
    imp.save(path)
    loaded = EMImputer.load(path)

    assert loaded._backend_name == "tensorflow"
    out_after = loaded.transform(X)
    np.testing.assert_allclose(out_before, out_after, rtol=1e-4, atol=1e-5)
