# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.ai.inversion — EMInverter1D, EnsembleInverter."""
import os
import tempfile
import unittest

import numpy as np


def _has_backend():
    from pycsamt.backends import get_backend
    return get_backend() != "none"


@unittest.skipUnless(_has_backend(), "no DL backend available")
class TestEMInverter1DWithBackend(unittest.TestCase):
    """Full fit/predict tests — skipped when no DL framework is installed."""

    @classmethod
    def setUpClass(cls):
        from pycsamt.forward.batch import generate_dataset
        cls.ds = generate_dataset(
            solver="MT1D",
            n_samples=120,
            freqs=np.logspace(0, 3, 20),
            n_layers=3,
            seed=42,
            verbose=False,
        )

    def test_fit_cnn1d(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D(arch="cnn1d", n_layers=3, solver="mt1d")
        inv.fit(self.ds, epochs=3, batch_size=32, verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_fit_resnet(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D(arch="resnet", n_layers=3, solver="mt1d")
        inv.fit(self.ds, epochs=3, batch_size=32, verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_predict_shape(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D(arch="cnn1d", n_layers=3)
        inv.fit(self.ds, epochs=3, batch_size=32, verbose=False)
        y_pred = inv.predict(self.ds.X[:10])
        self.assertEqual(y_pred.shape[0], 10)

    def test_predict_output_dtype(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D(arch="cnn1d", n_layers=3)
        inv.fit(self.ds, epochs=3, verbose=False)
        y_pred = inv.predict(self.ds.X[:5])
        self.assertEqual(y_pred.dtype.kind, "f")  # float

    def test_score_returns_float(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D(arch="cnn1d", n_layers=3)
        inv.fit(self.ds, epochs=3, verbose=False)
        score = inv.score(self.ds.X, self.ds.y)
        self.assertIsInstance(float(score), float)

    def test_save_load(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D(arch="cnn1d", n_layers=3)
        inv.fit(self.ds, epochs=3, verbose=False)
        y1 = inv.predict(self.ds.X[:5])

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "inv1d.npz")
            inv.save(path)
            inv2 = EMInverter1D.load(path)

        y2 = inv2.predict(self.ds.X[:5])
        np.testing.assert_allclose(y1, y2, rtol=1e-4)

    def test_predict_models_length(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D(arch="cnn1d", n_layers=3)
        inv.fit(self.ds, epochs=3, verbose=False)
        models = inv.predict_models(self.ds.X[:4])
        self.assertEqual(len(models), 4)

    def test_predict_before_fit_raises(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D(arch="cnn1d", n_layers=3)
        with self.assertRaises(RuntimeError):
            inv.predict(self.ds.X[:2])

    def test_history_populated(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D(arch="cnn1d", n_layers=3)
        inv.fit(self.ds, epochs=4, verbose=False)
        self.assertIn("train_loss", inv._history)
        self.assertGreater(len(inv._history["train_loss"]), 0)

    def test_ensemble_predict_shape(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        from pycsamt.ai.inversion import EnsembleInverter
        base = EMInverter1D(arch="cnn1d", n_layers=3)
        ens = EnsembleInverter(base_estimator=base, n_estimators=2)
        ens.fit(self.ds, epochs=3, verbose=False)
        mean, std = ens.predict_with_uncertainty(self.ds.X[:6])
        self.assertEqual(mean.shape[0], 6)
        self.assertEqual(std.shape[0], 6)


class TestEMInverter1DInterface(unittest.TestCase):
    """Interface and validation tests — no DL required."""

    def test_importable(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D()
        self.assertIsNotNone(inv)

    def test_default_arch(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D()
        self.assertEqual(inv.arch, "resnet")

    def test_repr_not_empty(self):
        from pycsamt.ai.inversion.inv1d import EMInverter1D
        inv = EMInverter1D(arch="cnn1d", n_layers=5)
        self.assertIn("EMInverter1D", repr(inv))

    def test_ai_package_imports(self):
        """Top-level pycsamt.ai should export all public classes."""
        from pycsamt.ai import (
            EMInverter1D, EnsembleInverter,
            CNN1DNet, ResNet1DNet, FCN1DNet,
            Normalizer, EMDataset, EMTrainer,
            EMDenoiser, AnomalyDetector, DimensionalityClassifier,
            EMQCScorer,
        )

    def test_list_pretrained_returns_list(self):
        from pycsamt.ai import list_pretrained
        zoo = list_pretrained()
        self.assertIsInstance(zoo, (list, dict))


class TestEMCheckpoint(unittest.TestCase):

    def test_save_load_roundtrip(self):
        from pycsamt.ai._base import EMCheckpoint
        params = {"arch": "cnn1d", "n_layers": 3}
        weights = {
            "w0": np.array([[1.0, 2.0], [3.0, 4.0]]),
            "b0": np.array([0.1, 0.2]),
        }
        history = {"train_loss": [0.5, 0.4, 0.3]}

        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "ckpt.npz")
            ckpt = EMCheckpoint(params=params, weights=weights, history=history)
            ckpt.save(path)
            ckpt2 = EMCheckpoint.load(path)

        self.assertEqual(ckpt2.params["arch"], "cnn1d")
        np.testing.assert_array_equal(ckpt2.weights["w0"], weights["w0"])
        self.assertEqual(ckpt2.history["train_loss"], history["train_loss"])

    def test_missing_file_raises(self):
        from pycsamt.ai._base import EMCheckpoint
        with self.assertRaises(FileNotFoundError):
            EMCheckpoint.load("/nonexistent/path.npz")


if __name__ == "__main__":
    unittest.main()
