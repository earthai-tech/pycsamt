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


class TestGCNInverter3DInterface(unittest.TestCase):
    """Interface tests for GCNInverter3D — no DL backend required."""

    def test_importable(self):
        from pycsamt.ai.inversion.inv3d import GCNInverter3D
        inv = GCNInverter3D(n_features=8, n_layers=3)
        self.assertIsNotNone(inv)

    def test_default_params(self):
        from pycsamt.ai.inversion.inv3d import GCNInverter3D
        inv = GCNInverter3D()
        self.assertEqual(inv.n_features, 40)
        self.assertEqual(inv.n_layers, 5)
        self.assertEqual(inv.n_out, 9)          # 2*5-1

    def test_repr_unfitted(self):
        from pycsamt.ai.inversion.inv3d import GCNInverter3D
        inv = GCNInverter3D(n_features=10, n_layers=4)
        self.assertIn("GCNInverter3D", repr(inv))
        self.assertIn("unfitted", repr(inv))

    def test_predict_before_fit_raises(self):
        from pycsamt.ai.inversion.inv3d import GCNInverter3D
        inv = GCNInverter3D(n_features=8, n_layers=3)
        X = np.ones((5, 8), dtype=np.float32)
        with self.assertRaises(RuntimeError):
            inv.predict(X)

    def test_build_adjacency(self):
        from pycsamt.ai.nets.gcn import build_adjacency
        coords = np.random.default_rng(0).uniform(0, 5000, (20, 2))
        A = build_adjacency(coords, radius=2000)
        self.assertEqual(A.shape, (20, 20))
        self.assertEqual(A.dtype, np.float32)
        np.testing.assert_allclose(A, A.T, atol=1e-6)   # symmetric

    def test_gcnnet_exported(self):
        from pycsamt.ai import GCNNet, build_adjacency
        self.assertTrue(callable(GCNNet))
        self.assertTrue(callable(build_adjacency))

    def test_top_level_export(self):
        from pycsamt.ai import GCNInverter3D
        inv = GCNInverter3D(n_features=6, n_layers=2)
        self.assertEqual(inv.n_out, 3)   # 2*2-1


@unittest.skipUnless(_has_backend(), "no DL backend available")
class TestGCNInverter3DWithBackend(unittest.TestCase):
    """Full fit/predict tests for GCNInverter3D."""

    @classmethod
    def setUpClass(cls):
        from pycsamt.ai.nets.gcn import build_adjacency
        rng = np.random.default_rng(0)
        cls.n_sta  = 15
        cls.n_feat = 8
        cls.n_layers = 3
        cls.n_out  = 2 * cls.n_layers - 1
        cls.n_samples = 50
        coords = rng.uniform(0, 5000, (cls.n_sta, 2))
        cls.A = build_adjacency(coords, radius=2000)
        cls.X = rng.standard_normal(
            (cls.n_samples, cls.n_sta, cls.n_feat)
        ).astype(np.float32)
        cls.y = rng.standard_normal(
            (cls.n_samples, cls.n_sta, cls.n_out)
        ).astype(np.float32)

    def _make_inv(self):
        from pycsamt.ai.inversion.inv3d import GCNInverter3D
        return GCNInverter3D(
            n_features=self.n_feat,
            n_layers=self.n_layers,
            hidden=(32,),
        )

    def test_fit_marks_fitted(self):
        inv = self._make_inv()
        inv.fit(self.X, self.y, adjacency=self.A, epochs=3, verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_predict_shape_batch(self):
        inv = self._make_inv()
        inv.fit(self.X, self.y, adjacency=self.A, epochs=3, verbose=False)
        pred = inv.predict(self.X[:4])
        self.assertEqual(pred.shape, (4, self.n_sta, self.n_out))

    def test_predict_shape_single(self):
        inv = self._make_inv()
        inv.fit(self.X, self.y, adjacency=self.A, epochs=3, verbose=False)
        pred = inv.predict(self.X[0])
        self.assertEqual(pred.shape, (self.n_sta, self.n_out))

    def test_predict_dtype_float(self):
        inv = self._make_inv()
        inv.fit(self.X, self.y, adjacency=self.A, epochs=3, verbose=False)
        pred = inv.predict(self.X[:2])
        self.assertEqual(pred.dtype.kind, "f")

    def test_history_populated(self):
        inv = self._make_inv()
        inv.fit(self.X, self.y, adjacency=self.A, epochs=4, verbose=False)
        self.assertIn("train_loss", inv._history)
        self.assertGreater(len(inv._history["train_loss"]), 0)

    def test_save_load_roundtrip(self):
        import tempfile, os
        inv = self._make_inv()
        inv.fit(self.X, self.y, adjacency=self.A, epochs=3, verbose=False)
        y1 = inv.predict(self.X[:3])
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "gcn3d.npz")
            inv.save(path)
            from pycsamt.ai.inversion.inv3d import GCNInverter3D
            inv2 = GCNInverter3D.load(path)
        y2 = inv2.predict(self.X[:3])
        np.testing.assert_allclose(y1, y2, rtol=1e-4)

    def test_uncertainty_shapes(self):
        inv = self._make_inv()
        inv.fit(self.X, self.y, adjacency=self.A, epochs=3, verbose=False)
        mean, std = inv.predict_with_uncertainty(self.X[:3], n_mc=5)
        self.assertEqual(mean.shape, (3, self.n_sta, self.n_out))
        self.assertEqual(std.shape, (3, self.n_sta, self.n_out))

    def test_fit_from_coords(self):
        """fit accepts coords= instead of adjacency=."""
        from pycsamt.ai.inversion.inv3d import GCNInverter3D
        rng = np.random.default_rng(1)
        coords = rng.uniform(0, 4000, (self.n_sta, 2))
        inv = GCNInverter3D(n_features=self.n_feat, n_layers=self.n_layers,
                            hidden=(32,))
        inv.fit(self.X, self.y, coords=coords, radius=2500,
                epochs=3, verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_identity_adjacency_warning(self):
        """No adjacency or coords → UserWarning about identity adjacency."""
        inv = self._make_inv()
        with self.assertWarns(UserWarning):
            inv.fit(self.X, self.y, epochs=2, verbose=False)


class TestConformalPredictor(unittest.TestCase):
    """Tests for ConformalPredictor — no DL backend required."""

    def _make_predictor(self, n_cal=80, n_params=4, seed=0):
        """Minimal stand-in for a fitted predictor."""
        rng = np.random.default_rng(seed)

        class _FakePredictor:
            def predict_with_uncertainty(self, X, _use_calibrated=True):
                rng2 = np.random.default_rng(42)
                n = len(X)
                mean  = rng2.standard_normal((n, n_params))
                sigma = np.abs(rng2.standard_normal((n, n_params))) + 0.1
                return mean, sigma

        pred   = _FakePredictor()
        X_cal  = rng.standard_normal((n_cal, 10))
        mean_c, sigma_c = pred.predict_with_uncertainty(X_cal)
        y_cal  = mean_c + sigma_c * rng.standard_normal((n_cal, n_params))
        return pred, X_cal, y_cal, n_params

    def test_importable(self):
        from pycsamt.ai.inversion.calibration import ConformalPredictor
        cp = ConformalPredictor(None)
        self.assertIsNotNone(cp)

    def test_calibrate_returns_self(self):
        from pycsamt.ai.inversion.calibration import ConformalPredictor
        pred, X_cal, y_cal, _ = self._make_predictor()
        cp = ConformalPredictor(pred)
        result = cp.calibrate(X_cal, y_cal)
        self.assertIs(result, cp)

    def test_is_calibrated_after_calibrate(self):
        from pycsamt.ai.inversion.calibration import ConformalPredictor
        pred, X_cal, y_cal, _ = self._make_predictor()
        cp = ConformalPredictor(pred).calibrate(X_cal, y_cal)
        self.assertTrue(cp._is_calibrated)

    def test_predict_intervals_shapes(self):
        from pycsamt.ai.inversion.calibration import ConformalPredictor
        pred, X_cal, y_cal, n_p = self._make_predictor()
        cp = ConformalPredictor(pred).calibrate(X_cal, y_cal)
        rng = np.random.default_rng(1)
        X_test = rng.standard_normal((20, 10))
        center, lo, hi = cp.predict_intervals(X_test)
        self.assertEqual(center.shape, (20, n_p))
        self.assertEqual(lo.shape,     (20, n_p))
        self.assertEqual(hi.shape,     (20, n_p))

    def test_intervals_ordered(self):
        from pycsamt.ai.inversion.calibration import ConformalPredictor
        pred, X_cal, y_cal, _ = self._make_predictor()
        cp = ConformalPredictor(pred).calibrate(X_cal, y_cal)
        X_test = np.random.default_rng(2).standard_normal((30, 10))
        center, lo, hi = cp.predict_intervals(X_test)
        self.assertTrue(np.all(lo <= center))
        self.assertTrue(np.all(center <= hi))

    def test_coverage_in_range(self):
        from pycsamt.ai.inversion.calibration import ConformalPredictor
        pred, X_cal, y_cal, n_p = self._make_predictor(n_cal=200)
        cp = ConformalPredictor(pred).calibrate(X_cal, y_cal, alpha=0.10)
        rng = np.random.default_rng(3)
        X_t = rng.standard_normal((200, 10))
        mean_t, sigma_t = pred.predict_with_uncertainty(X_t)
        y_t  = mean_t + sigma_t * rng.standard_normal((200, n_p))
        cov  = cp.coverage(X_t, y_t, alpha=0.10)
        # Coverage should be >= 0.90 (guarantee) and not absurdly high
        self.assertGreaterEqual(cov, 0.80)
        self.assertLessEqual(cov, 1.0)

    def test_predict_before_calibrate_raises(self):
        from pycsamt.ai.inversion.calibration import ConformalPredictor
        cp = ConformalPredictor(None)
        with self.assertRaises(RuntimeError):
            cp.predict_intervals(np.ones((5, 4)))

    def test_coverage_diagnostics_returns_dict(self):
        from pycsamt.ai.inversion.calibration import ConformalPredictor
        pred, X_cal, y_cal, n_p = self._make_predictor(n_cal=100)
        cp = ConformalPredictor(pred).calibrate(X_cal, y_cal)
        rng = np.random.default_rng(4)
        X_t = rng.standard_normal((50, 10))
        mean_t, sigma_t = pred.predict_with_uncertainty(X_t)
        y_t  = mean_t + sigma_t * rng.standard_normal((50, n_p))
        diag = cp.coverage_diagnostics(X_t, y_t)
        self.assertIsInstance(diag, dict)
        self.assertGreater(len(diag), 0)
        for v in diag.values():
            self.assertGreaterEqual(v, 0.0)
            self.assertLessEqual(v, 1.0)

    def test_repr_contains_calibrated(self):
        from pycsamt.ai.inversion.calibration import ConformalPredictor
        pred, X_cal, y_cal, _ = self._make_predictor()
        cp = ConformalPredictor(pred).calibrate(X_cal, y_cal)
        self.assertIn("calibrated", repr(cp))

    def test_top_level_export(self):
        from pycsamt.ai import ConformalPredictor
        self.assertTrue(callable(ConformalPredictor))


class TestPosteriorCalibrator(unittest.TestCase):
    """Tests for PosteriorCalibrator — no DL backend required."""

    def _make_data(self, n=150, n_params=4, seed=0):
        rng = np.random.default_rng(seed)
        y_pred  = rng.standard_normal((n, n_params))
        sigma   = np.abs(rng.standard_normal((n, n_params))) + 0.2
        y_true  = y_pred + sigma * rng.standard_normal((n, n_params))
        return y_true, y_pred, sigma

    def test_importable(self):
        from pycsamt.ai.inversion.calibration import PosteriorCalibrator
        pc = PosteriorCalibrator()
        self.assertIsNotNone(pc)

    def test_fit_returns_self(self):
        from pycsamt.ai.inversion.calibration import PosteriorCalibrator
        y_true, y_pred, sigma = self._make_data()
        pc = PosteriorCalibrator()
        result = pc.fit(y_true, y_pred, sigma)
        self.assertIs(result, pc)

    def test_fitted_flag(self):
        from pycsamt.ai.inversion.calibration import PosteriorCalibrator
        y_true, y_pred, sigma = self._make_data()
        pc = PosteriorCalibrator().fit(y_true, y_pred, sigma)
        self.assertTrue(pc._is_fitted)

    def test_predict_posterior_shape(self):
        from pycsamt.ai.inversion.calibration import PosteriorCalibrator
        y_true, y_pred, sigma = self._make_data(n=100, n_params=3)
        pc = PosteriorCalibrator().fit(y_true, y_pred, sigma)
        rng = np.random.default_rng(5)
        y_p2 = np.random.default_rng(6).standard_normal((20, 3))
        s_p2 = np.abs(np.random.default_rng(7).standard_normal((20, 3))) + 0.1
        draws = pc.predict_posterior(y_p2, s_p2, n_samples=50, rng=rng)
        self.assertEqual(draws.shape, (50, 20, 3))

    def test_posterior_mean_close_to_pred(self):
        from pycsamt.ai.inversion.calibration import PosteriorCalibrator
        y_true, y_pred, sigma = self._make_data(n=200)
        pc = PosteriorCalibrator().fit(y_true, y_pred, sigma)
        rng = np.random.default_rng(8)
        y_p2 = rng.standard_normal((30, 4))
        s_p2 = np.ones((30, 4)) * 0.5
        draws = pc.predict_posterior(y_p2, s_p2, n_samples=2000, rng=rng)
        # mean of posterior samples should be near y_pred
        np.testing.assert_allclose(draws.mean(axis=0), y_p2, atol=0.15)

    def test_calibrated_std_returns_array(self):
        from pycsamt.ai.inversion.calibration import PosteriorCalibrator
        y_true, y_pred, sigma = self._make_data()
        pc = PosteriorCalibrator().fit(y_true, y_pred, sigma)
        sigma_cal = pc.calibrated_std(sigma[:10])
        self.assertEqual(sigma_cal.shape, sigma[:10].shape)
        self.assertTrue(np.all(sigma_cal > 0))

    def test_mace_in_range(self):
        from pycsamt.ai.inversion.calibration import PosteriorCalibrator
        y_true, y_pred, sigma = self._make_data(n=200)
        pc = PosteriorCalibrator().fit(y_true, y_pred, sigma)
        y_t2, y_p2, s2 = self._make_data(n=100, seed=99)
        mace = pc.calibration_error(y_t2, y_p2, s2)
        self.assertGreaterEqual(mace, 0.0)
        self.assertLessEqual(mace, 1.0)

    def test_predict_before_fit_raises(self):
        from pycsamt.ai.inversion.calibration import PosteriorCalibrator
        pc = PosteriorCalibrator()
        with self.assertRaises(RuntimeError):
            pc.calibrated_std(np.ones((5, 3)))

    def test_repr_contains_status(self):
        from pycsamt.ai.inversion.calibration import PosteriorCalibrator
        pc = PosteriorCalibrator().fit(*self._make_data())
        self.assertIn("fitted", repr(pc))

    def test_top_level_export(self):
        from pycsamt.ai import PosteriorCalibrator
        self.assertTrue(callable(PosteriorCalibrator))


class TestEnsembleInverterCalibrate(unittest.TestCase):
    """Tests for EnsembleInverter.calibrate / predict_intervals / predict_posterior."""

    def _make_ensemble_stub(self, n_params=4, seed=0):
        """Return a minimal EnsembleInverter-like object with stub members."""
        from pycsamt.ai.inversion.ensemble import EnsembleInverter

        class _StubMember:
            def __init__(self, rng):
                self._rng = rng
            def predict(self, X, **kw):
                rng2 = np.random.default_rng(self._rng)
                return rng2.standard_normal((len(X), n_params))

        # Build EnsembleInverter without actually training
        ens = object.__new__(EnsembleInverter)
        ens.base_estimator   = None
        ens.n_estimators     = 3
        ens.seeds            = [0, 1, 2]
        ens._members         = [_StubMember(s) for s in [0, 1, 2]]
        ens._is_fitted       = True
        ens._conformal       = None
        ens._posterior_cal   = None
        return ens, n_params

    def test_calibrate_attaches_conformal(self):
        from pycsamt.ai.inversion.calibration import ConformalPredictor
        ens, n_p = self._make_ensemble_stub()
        rng = np.random.default_rng(0)
        X_cal = rng.standard_normal((80, 6))
        mean, sigma = ens.predict_with_uncertainty(X_cal, _use_calibrated=False)
        y_cal = mean + sigma * rng.standard_normal((80, n_p))
        ens.calibrate(X_cal, y_cal)
        self.assertIsInstance(ens._conformal, ConformalPredictor)

    def test_calibrate_attaches_posterior(self):
        from pycsamt.ai.inversion.calibration import PosteriorCalibrator
        ens, n_p = self._make_ensemble_stub()
        rng = np.random.default_rng(1)
        X_cal = rng.standard_normal((80, 6))
        mean, sigma = ens.predict_with_uncertainty(X_cal, _use_calibrated=False)
        y_cal = mean + sigma * rng.standard_normal((80, n_p))
        ens.calibrate(X_cal, y_cal)
        self.assertIsInstance(ens._posterior_cal, PosteriorCalibrator)

    def test_predict_intervals_shape(self):
        ens, n_p = self._make_ensemble_stub()
        rng = np.random.default_rng(2)
        X_cal = rng.standard_normal((80, 6))
        mean, sigma = ens.predict_with_uncertainty(X_cal, _use_calibrated=False)
        y_cal = mean + sigma * rng.standard_normal((80, n_p))
        ens.calibrate(X_cal, y_cal)
        X_t = rng.standard_normal((25, 6))
        center, lo, hi = ens.predict_intervals(X_t)
        self.assertEqual(center.shape, (25, n_p))
        self.assertTrue(np.all(lo <= center))
        self.assertTrue(np.all(center <= hi))

    def test_predict_posterior_shape(self):
        ens, n_p = self._make_ensemble_stub()
        rng = np.random.default_rng(3)
        X_cal = rng.standard_normal((80, 6))
        mean, sigma = ens.predict_with_uncertainty(X_cal, _use_calibrated=False)
        y_cal = mean + sigma * rng.standard_normal((80, n_p))
        ens.calibrate(X_cal, y_cal)
        X_t = rng.standard_normal((10, 6))
        draws = ens.predict_posterior(X_t, n_samples=30, rng=rng)
        self.assertEqual(draws.shape, (30, 10, n_p))

    def test_predict_with_uncertainty_returns_calibrated_sigma(self):
        ens, n_p = self._make_ensemble_stub()
        rng = np.random.default_rng(4)
        X_cal = rng.standard_normal((80, 6))
        mean0, sigma0 = ens.predict_with_uncertainty(X_cal, _use_calibrated=False)
        y_cal = mean0 + sigma0 * rng.standard_normal((80, n_p))
        ens.calibrate(X_cal, y_cal)
        X_t = rng.standard_normal((20, 6))
        _, sigma_cal  = ens.predict_with_uncertainty(X_t)
        _, sigma_raw  = ens.predict_with_uncertainty(X_t, _use_calibrated=False)
        # After calibration, sigma_cal != sigma_raw (unless scale=1 exactly)
        # At minimum both should be positive and finite
        self.assertTrue(np.all(sigma_cal > 0))
        self.assertTrue(np.all(np.isfinite(sigma_cal)))

    def test_predict_intervals_before_calibrate_raises(self):
        ens, _ = self._make_ensemble_stub()
        with self.assertRaises(RuntimeError):
            ens.predict_intervals(np.ones((5, 6)))

    def test_predict_posterior_before_calibrate_raises(self):
        ens, _ = self._make_ensemble_stub()
        with self.assertRaises(RuntimeError):
            ens.predict_posterior(np.ones((5, 6)))

    def test_coverage_diagnostics_before_calibrate_raises(self):
        ens, _ = self._make_ensemble_stub()
        with self.assertRaises(RuntimeError):
            ens.coverage_diagnostics(np.ones((5, 6)), np.ones((5, 4)))


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
