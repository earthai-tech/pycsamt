# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.ai.processing — EMDenoiser, AnomalyDetector, DimensionalityClassifier."""

import unittest

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# EMDenoiser
# ─────────────────────────────────────────────────────────────────────────────


class TestEMDenoiser(unittest.TestCase):
    def _make_X(self, n=60, n_comp=4, n_freqs=24):
        rng = np.random.default_rng(7)
        return rng.standard_normal((n, n_comp, n_freqs)).astype(np.float32)

    def test_no_args_init(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        den = EMDenoiser()
        self.assertIsNone(den.n_freqs)

    def test_infers_n_freqs(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        X = self._make_X(n_freqs=32)
        den = EMDenoiser()
        den.fit(X, epochs=2, verbose=False)
        self.assertEqual(den.n_freqs, 32)

    def test_explicit_n_freqs(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        X = self._make_X(n_freqs=16)
        den = EMDenoiser(n_freqs=16)
        den.fit(X, epochs=2, verbose=False)
        self.assertEqual(den.n_freqs, 16)

    def test_n_freqs_mismatch_raises(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        X = self._make_X(n_freqs=24)
        den = EMDenoiser(n_freqs=32)
        with self.assertRaises(ValueError):
            den.fit(X, epochs=1, verbose=False)

    def test_transform_shape(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        X = self._make_X()
        den = EMDenoiser()
        den.fit(X, epochs=2, verbose=False)
        X2 = den.transform(X)
        self.assertEqual(X2.shape, X.shape)

    def test_transform_before_fit_raises(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        X = self._make_X()
        den = EMDenoiser()
        with self.assertRaises(RuntimeError):
            den.transform(X)

    def test_input_must_be_3d(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        den = EMDenoiser()
        with self.assertRaises(ValueError):
            den.fit(np.ones((10, 4)), epochs=1, verbose=False)

    def test_fit_transform_dtype(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        X = self._make_X()
        den = EMDenoiser()
        den.fit(X, epochs=2, verbose=False)
        X2 = den.transform(X)
        self.assertEqual(X2.dtype, np.float32)

    def test_output_finite(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        X = self._make_X()
        den = EMDenoiser()
        den.fit(X, epochs=2, verbose=False)
        X2 = den.transform(X)
        self.assertTrue(np.all(np.isfinite(X2)))

    def test_repr_before_fit(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        den = EMDenoiser()
        self.assertIn("unfitted", repr(den))

    def test_repr_after_fit(self):
        from pycsamt.ai.processing.denoise import EMDenoiser

        X = self._make_X()
        den = EMDenoiser()
        den.fit(X, epochs=1, verbose=False)
        self.assertIn("fitted", repr(den))


# ─────────────────────────────────────────────────────────────────────────────
# AnomalyDetector
# ─────────────────────────────────────────────────────────────────────────────


class TestAnomalyDetector(unittest.TestCase):
    def _make_X(self, n=120, n_feat=40):
        rng = np.random.default_rng(3)
        return rng.standard_normal((n, n_feat)).astype(np.float32)

    def test_no_args_init(self):
        from pycsamt.ai.processing.anomaly import (
            AnomalyDetector,
        )

        det = AnomalyDetector()
        self.assertIsNone(det.n_features)

    def test_infers_n_features(self):
        from pycsamt.ai.processing.anomaly import (
            AnomalyDetector,
        )

        X = self._make_X(n_feat=48)
        det = AnomalyDetector()
        det.fit(X, epochs=3, verbose=False)
        self.assertEqual(det.n_features, 48)

    def test_explicit_n_features(self):
        from pycsamt.ai.processing.anomaly import (
            AnomalyDetector,
        )

        X = self._make_X(n_feat=40)
        det = AnomalyDetector(n_features=40)
        det.fit(X, epochs=3, verbose=False)
        self.assertEqual(det.n_features, 40)

    def test_n_features_mismatch_raises(self):
        from pycsamt.ai.processing.anomaly import (
            AnomalyDetector,
        )

        X = self._make_X(n_feat=40)
        det = AnomalyDetector(n_features=50)
        with self.assertRaises(ValueError):
            det.fit(X, epochs=2, verbose=False)

    def test_transform_shape(self):
        from pycsamt.ai.processing.anomaly import (
            AnomalyDetector,
        )

        X = self._make_X()
        det = AnomalyDetector()
        det.fit(X, epochs=3, verbose=False)
        scores = det.transform(X)
        self.assertEqual(scores.shape, (len(X),))

    def test_scores_non_negative(self):
        from pycsamt.ai.processing.anomaly import (
            AnomalyDetector,
        )

        X = self._make_X()
        det = AnomalyDetector()
        det.fit(X, epochs=3, verbose=False)
        scores = det.transform(X)
        self.assertTrue(np.all(scores >= 0))

    def test_flag_anomalies_shape(self):
        from pycsamt.ai.processing.anomaly import (
            AnomalyDetector,
        )

        X = self._make_X()
        det = AnomalyDetector()
        det.fit(X, epochs=3, verbose=False)
        flags = det.flag_anomalies(X)
        self.assertEqual(flags.shape, (len(X),))
        self.assertEqual(flags.dtype, bool)

    def test_flag_percentile(self):
        """5% of normal training data should be flagged at 95th percentile."""
        from pycsamt.ai.processing.anomaly import (
            AnomalyDetector,
        )

        X = self._make_X(n=200)
        det = AnomalyDetector(threshold_percentile=95.0)
        det.fit(X, epochs=5, verbose=False)
        flags = det.flag_anomalies(X)
        frac = flags.mean()
        # Approximately 5% flagged; allow wide tolerance for small datasets
        self.assertLess(frac, 0.20)

    def test_transform_before_fit_raises(self):
        from pycsamt.ai.processing.anomaly import (
            AnomalyDetector,
        )

        X = self._make_X()
        det = AnomalyDetector()
        with self.assertRaises(RuntimeError):
            det.transform(X)

    def test_nan_in_X(self):
        from pycsamt.ai.processing.anomaly import (
            AnomalyDetector,
        )

        X = self._make_X()
        X[0, 0] = np.nan
        det = AnomalyDetector()
        det.fit(X, epochs=3, verbose=False)
        scores = det.transform(X)
        self.assertTrue(np.all(np.isfinite(scores)))


# ─────────────────────────────────────────────────────────────────────────────
# DimensionalityClassifier
# ─────────────────────────────────────────────────────────────────────────────


class TestDimensionalityClassifier(unittest.TestCase):
    def _make_X(self, n=300):
        rng = np.random.default_rng(11)
        # Feature vector: [beta_abs, ellipt_abs, logrho_det, phi_det, tip_amp]
        X = rng.standard_normal((n, 5)).astype(np.float32)
        X[:, 0] = np.abs(X[:, 0]) * 4  # beta_abs ≥ 0
        X[:, 1] = np.abs(X[:, 1]) * 0.3  # ellipt_abs ≥ 0
        return X

    def test_fit_no_y(self):
        from pycsamt.ai.processing.classify import (
            DimensionalityClassifier,
        )

        X = self._make_X()
        clf = DimensionalityClassifier()
        clf.fit(X, verbose=False)
        self.assertTrue(clf._is_fitted)

    def test_fit_with_y(self):
        from pycsamt.ai.processing.classify import (
            DimensionalityClassifier,
        )

        X = self._make_X()
        y = np.random.randint(0, 3, len(X))
        clf = DimensionalityClassifier()
        clf.fit(X, y, verbose=False)
        self.assertTrue(clf._is_fitted)

    def test_predict_shape(self):
        from pycsamt.ai.processing.classify import (
            DimensionalityClassifier,
        )

        X = self._make_X()
        clf = DimensionalityClassifier()
        clf.fit(X, verbose=False)
        preds = clf.predict(X)
        self.assertEqual(preds.shape, (len(X),))

    def test_predict_labels_valid(self):
        from pycsamt.ai.processing.classify import (
            DimensionalityClassifier,
        )

        X = self._make_X()
        clf = DimensionalityClassifier()
        clf.fit(X, verbose=False)
        preds = clf.predict(X)
        self.assertTrue(np.all(np.isin(preds, [0, 1, 2])))

    def test_transform_proba_shape(self):
        from pycsamt.ai.processing.classify import (
            DimensionalityClassifier,
        )

        X = self._make_X()
        clf = DimensionalityClassifier()
        clf.fit(X, verbose=False)
        proba = clf.transform(X)
        self.assertEqual(proba.shape, (len(X), 3))

    def test_proba_sums_to_one(self):
        from pycsamt.ai.processing.classify import (
            DimensionalityClassifier,
        )

        X = self._make_X()
        clf = DimensionalityClassifier()
        clf.fit(X, verbose=False)
        proba = clf.transform(X)
        np.testing.assert_allclose(proba.sum(axis=1), 1.0, atol=1e-5)

    def test_predict_strike_shape(self):
        from pycsamt.ai.processing.classify import (
            DimensionalityClassifier,
        )

        X = self._make_X()
        clf = DimensionalityClassifier()
        clf.fit(X, verbose=False)
        strike = clf.predict_strike(X)
        self.assertEqual(strike.shape, (len(X),))

    def test_transform_before_fit_raises(self):
        from pycsamt.ai.processing.classify import (
            DimensionalityClassifier,
        )

        X = self._make_X()
        clf = DimensionalityClassifier()
        with self.assertRaises(RuntimeError):
            clf.transform(X)

    def test_rule_labels_consistency(self):
        """Rule-based labels should agree with the heuristic thresholds."""
        from pycsamt.ai.processing.classify import (
            _rule_labels,
        )

        beta = np.array([1.0, 5.0, 2.0, 8.0])
        ellipt = np.array([0.1, 0.1, 0.5, 0.5])
        labels = _rule_labels(beta, ellipt, skew_th=3.0, ellipt_th=0.2)
        np.testing.assert_array_equal(labels, [0, 2, 1, 2])

    def test_dataframe_input(self):
        import pandas as pd

        from pycsamt.ai.processing.classify import (
            DimensionalityClassifier,
        )

        X = self._make_X(n=100)
        df = pd.DataFrame(
            X,
            columns=[
                "beta_abs",
                "ellipt_abs",
                "logrho_det",
                "phi_det",
                "tip_amp",
            ],
        )
        clf = DimensionalityClassifier()
        clf.fit(df, verbose=False)
        preds = clf.predict(df)
        self.assertEqual(preds.shape, (100,))


# ─────────────────────────────────────────────────────────────────────────────
# EMQCScorer  (interface smoke test — no emtools dependency)
# ─────────────────────────────────────────────────────────────────────────────


class TestEMQCScorerImport(unittest.TestCase):
    def test_importable(self):
        from pycsamt.ai.processing.qc import EMQCScorer

        obj = EMQCScorer()
        self.assertIsNotNone(obj)

    def test_has_fit_and_transform(self):
        from pycsamt.ai.processing.qc import EMQCScorer

        obj = EMQCScorer()
        self.assertTrue(hasattr(obj, "fit"))
        self.assertTrue(hasattr(obj, "transform"))


if __name__ == "__main__":
    unittest.main()
