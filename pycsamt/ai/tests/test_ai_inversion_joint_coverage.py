# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Coverage tests for pycsamt.ai.inversion.joint.JointInverter.

Mirrors the fixture/test-writing conventions of
``test_inversion.py`` (EMInverter1D, GCNInverter3D): real backend
training with small epoch counts, synthetic multi-modality feature
arrays built like :func:`pycsamt.forward.batch.generate_dataset`
output, and save/load round trips.
"""

from __future__ import annotations

import unittest

import numpy as np


def _has_backend():
    from pycsamt.backends import get_backend

    return get_backend() != "none"


def _synthetic_joint_data(
    n=60, n_features_list=(12, 8), n_layers=3, seed=0
):
    rng = np.random.default_rng(seed)
    X_list = [
        rng.standard_normal((n, nf)).astype(np.float32)
        for nf in n_features_list
    ]
    n_out = 2 * n_layers - 1
    y = rng.uniform(1.0, 3.0, size=(n, n_out)).astype(np.float32)
    # thickness columns must be positive for the log10 transform
    y[:, n_layers:] = rng.uniform(10.0, 200.0, size=(n, n_layers - 1))
    return X_list, y


@unittest.skipUnless(_has_backend(), "no DL backend available")
class TestJointInverterConstruction(unittest.TestCase):
    def test_defaults(self):
        from pycsamt.ai.inversion.joint import JointInverter

        inv = JointInverter(n_features_list=(120, 48), n_layers=5)
        self.assertEqual(inv.n_features_list, (120, 48))
        self.assertEqual(inv.n_layers, 5)
        self.assertEqual(inv._n_out, 9)
        self.assertFalse(inv._is_fitted)

    def test_repr_unfitted(self):
        from pycsamt.ai.inversion.joint import JointInverter

        inv = JointInverter(n_features_list=(10, 10))
        self.assertIn("unfitted", repr(inv))
        self.assertIn("modalities=2", repr(inv))

    def test_predict_before_fit_raises(self):
        from pycsamt.ai.inversion.joint import JointInverter

        inv = JointInverter(n_features_list=(10, 10))
        with self.assertRaises(RuntimeError):
            inv.predict([np.zeros((2, 10)), np.zeros((2, 10))])


@unittest.skipUnless(_has_backend(), "no DL backend available")
class TestJointInverterFitPredictTorch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X_list, cls.y = _synthetic_joint_data()

    def _new(self, **kwargs):
        from pycsamt.ai.inversion.joint import JointInverter

        return JointInverter(
            n_features_list=(12, 8),
            n_layers=3,
            hidden_dim=16,
            growth_rate=8,
            n_dense_layers=2,
            **kwargs,
        )

    def test_fit_marks_fitted_and_history(self):
        inv = self._new()
        inv.fit(self.X_list, self.y, epochs=3, batch_size=16, verbose=False)
        self.assertTrue(inv._is_fitted)
        self.assertIn("train_loss", inv._history)
        self.assertIn("val_loss", inv._history)
        self.assertIn("best_val_loss", inv._meta)
        self.assertEqual(inv._backend_name, "torch")
        self.assertIn("fitted", repr(inv))

    def test_predict_shape_log_rho_default(self):
        inv = self._new()
        inv.fit(self.X_list, self.y, epochs=3, batch_size=16, verbose=False)
        y_pred = inv.predict(self.X_list)
        self.assertEqual(y_pred.shape, (60, 5))
        self.assertTrue(np.all(np.isfinite(y_pred)))

    def test_predict_as_log_rho_false_gives_linear_resistivity(self):
        inv = self._new()
        inv.fit(self.X_list, self.y, epochs=2, batch_size=16, verbose=False)
        y_log = inv.predict(self.X_list, as_log_rho=True)
        y_lin = inv.predict(self.X_list, as_log_rho=False)
        np.testing.assert_allclose(
            10.0 ** y_log[:, :3], y_lin[:, :3], rtol=1e-4
        )
        np.testing.assert_allclose(y_log[:, 3:], y_lin[:, 3:], rtol=1e-4)

    def test_predict_single_ndarray_wrapped_in_list(self):
        from pycsamt.ai.inversion.joint import JointInverter

        rng = np.random.default_rng(1)
        X = rng.standard_normal((20, 10)).astype(np.float32)
        y = rng.uniform(1.0, 3.0, size=(20, 5)).astype(np.float32)
        y[:, 3:] = rng.uniform(10.0, 100.0, size=(20, 2))

        inv = JointInverter(
            n_features_list=(10,),
            n_layers=3,
            hidden_dim=16,
            growth_rate=8,
            n_dense_layers=2,
        )
        inv.fit([X], y, epochs=2, batch_size=8, verbose=False)
        y_pred_list = inv.predict([X])
        y_pred_bare = inv.predict(X)
        np.testing.assert_allclose(y_pred_list, y_pred_bare)

    def test_log_thickness_false_skips_log_transform(self):
        inv = self._new(log_thickness=False)
        inv.fit(self.X_list, self.y, epochs=2, batch_size=16, verbose=False)
        y_pred = inv.predict(self.X_list)
        self.assertEqual(y_pred.shape, (60, 5))
        self.assertTrue(np.all(np.isfinite(y_pred)))

    def test_verbose_progress_and_early_stopping_prints(self):
        inv = self._new()
        # lr=0 freezes the network -> val loss is bit-identical every
        # epoch after the first, guaranteeing no_improve reaches
        # patience=1 deterministically (real early-stopping branch).
        inv.fit(
            self.X_list,
            self.y,
            epochs=10,
            batch_size=16,
            lr=0.0,
            patience=1,
            verbose=True,
        )
        self.assertLess(len(inv._history["train_loss"]), 10)

    def test_grad_clip_none_skips_clipping_branch(self):
        inv = self._new()
        inv.fit(
            self.X_list,
            self.y,
            epochs=2,
            batch_size=16,
            grad_clip=None,
            verbose=False,
        )
        self.assertTrue(inv._is_fitted)

    def test_three_modalities(self):
        from pycsamt.ai.inversion.joint import JointInverter

        X_list, y = _synthetic_joint_data(n=40, n_features_list=(10, 6, 4))
        inv = JointInverter(
            n_features_list=(10, 6, 4),
            n_layers=3,
            hidden_dim=16,
            growth_rate=8,
            n_dense_layers=2,
        )
        inv.fit(X_list, y, epochs=2, batch_size=8, verbose=False)
        y_pred = inv.predict(X_list)
        self.assertEqual(y_pred.shape, (40, 5))

    def test_save_load_round_trip(self):
        import tempfile
        from pathlib import Path

        from pycsamt.ai.inversion.joint import JointInverter

        inv = self._new()
        inv.fit(self.X_list, self.y, epochs=2, batch_size=16, verbose=False)
        y_before = inv.predict(self.X_list)

        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "joint.npz"
            inv.save(path)
            loaded = JointInverter.load(path)

        self.assertEqual(loaded._backend_name, "torch")
        y_after = loaded.predict(self.X_list)
        np.testing.assert_allclose(y_before, y_after, rtol=1e-4, atol=1e-5)
        np.testing.assert_allclose(
            loaded._meta["best_val_loss"], inv._meta["best_val_loss"]
        )


if __name__ == "__main__":
    unittest.main()
