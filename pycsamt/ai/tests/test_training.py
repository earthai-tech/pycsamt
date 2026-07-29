# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.ai training utilities — Normalizer, EMDataset, augmenters, metrics."""

import unittest

import numpy as np


def _has_backend():
    from pycsamt.backends import get_backend

    return get_backend() != "none"


# ─────────────────────────────────────────────────────────────────────────────
# Normalizer
# ─────────────────────────────────────────────────────────────────────────────


class TestNormalizer(unittest.TestCase):
    def setUp(self):
        from pycsamt.ai.training.dataset import Normalizer

        self.rng = np.random.default_rng(0)
        self.X = self.rng.standard_normal((100, 12)).astype(np.float32)
        self.norm = Normalizer()

    def test_fit_transform(self):
        Xn = self.norm.fit(self.X).transform(self.X)
        np.testing.assert_allclose(Xn.mean(axis=0), 0.0, atol=1e-4)
        np.testing.assert_allclose(Xn.std(axis=0), 1.0, atol=1e-4)

    def test_fit_transform_shortcut(self):
        from pycsamt.ai.training.dataset import Normalizer

        norm = Normalizer()
        Xn = norm.fit_transform(self.X)
        np.testing.assert_allclose(Xn.mean(axis=0), 0.0, atol=1e-4)

    def test_inverse_transform(self):
        Xn = self.norm.fit(self.X).transform(self.X)
        X_back = self.norm.inverse_transform(Xn)
        np.testing.assert_allclose(X_back, self.X, atol=1e-5)

    def test_serialisation_roundtrip(self):
        from pycsamt.ai.training.dataset import Normalizer

        self.norm.fit(self.X)
        d = self.norm.to_dict()
        norm2 = Normalizer.from_dict(d)
        Xn1 = self.norm.transform(self.X)
        Xn2 = norm2.transform(self.X)
        np.testing.assert_allclose(Xn1, Xn2, atol=1e-6)

    def test_nan_handling(self):
        X_with_nan = self.X.copy()
        X_with_nan[0, 0] = np.nan
        Xn = self.norm.fit(X_with_nan).transform(X_with_nan)
        # Non-NaN entries should be finite
        finite_mask = ~np.isnan(X_with_nan)
        self.assertTrue(np.isfinite(Xn[finite_mask]).all())


# ─────────────────────────────────────────────────────────────────────────────
# EMDataset (requires forward data)
# ─────────────────────────────────────────────────────────────────────────────


class TestEMDataset(unittest.TestCase):
    def _make_forward_ds(self, n=30, n_freqs=10, n_layers=3):
        from pycsamt.forward.batch import generate_dataset

        return generate_dataset(
            solver="MT1D",
            n_samples=n,
            freqs=np.logspace(0, 3, n_freqs),
            n_layers=n_layers,
            seed=7,
            verbose=False,
        )

    def test_length(self):
        from pycsamt.ai.training.dataset import EMDataset

        fds = self._make_forward_ds()
        ds = EMDataset(fds)
        self.assertEqual(len(ds), len(fds.X))

    @unittest.skipUnless(
        __import__("importlib").util.find_spec("torch") is not None,
        "torch required for EMDataset.__getitem__",
    )
    def test_getitem_shapes(self):
        from pycsamt.ai.training.dataset import EMDataset

        fds = self._make_forward_ds(n_layers=3, n_freqs=10)
        ds = EMDataset(fds)
        x, y = ds[0]
        self.assertEqual(x.ndim, 1)
        self.assertEqual(y.ndim, 1)

    def test_split_sizes(self):
        from pycsamt.ai.training.dataset import EMDataset

        fds = self._make_forward_ds(n=40)
        ds = EMDataset(fds)
        trn, val = ds.split(val_frac=0.25, seed=0)
        self.assertEqual(len(trn) + len(val), 40)

    def test_normalizers_preserved_after_split(self):
        from pycsamt.ai.training.dataset import EMDataset

        fds = self._make_forward_ds(n=40)
        ds = EMDataset(fds)
        trn, val = ds.split(val_frac=0.25, seed=0)
        # After split, both subsets share the same normaliser objects
        self.assertIs(trn.x_norm, val.x_norm)
        self.assertIs(trn.y_norm, val.y_norm)

    def test_x_shape(self):
        from pycsamt.ai.training.dataset import EMDataset

        fds = self._make_forward_ds(n=20, n_freqs=8)
        ds = EMDataset(fds)
        self.assertEqual(ds.X.shape[0], 20)

    def test_y_shape(self):
        from pycsamt.ai.training.dataset import EMDataset

        fds = self._make_forward_ds(n=20, n_layers=3)
        ds = EMDataset(fds)
        # n_layers=3: 3 rho + 2 thick = 5 params
        self.assertEqual(ds.y.shape[1], 5)


# ─────────────────────────────────────────────────────────────────────────────
# Augmenters
# ─────────────────────────────────────────────────────────────────────────────


class TestAugmentNoise(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(42)
        self.X = self.rng.standard_normal((50, 24)).astype(np.float32)
        self.y = self.rng.standard_normal((50, 7)).astype(np.float32)

    def test_output_shape(self):
        from pycsamt.ai.training.augment import AugmentNoise

        aug = AugmentNoise(sigma=0.05)
        Xa, ya = aug(self.X, self.y)
        self.assertEqual(Xa.shape, self.X.shape)
        self.assertEqual(ya.shape, self.y.shape)

    def test_modifies_X(self):
        from pycsamt.ai.training.augment import AugmentNoise

        aug = AugmentNoise(sigma=0.5)
        Xa, _ = aug(self.X, self.y)
        self.assertFalse(np.allclose(Xa, self.X))

    def test_y_unchanged(self):
        from pycsamt.ai.training.augment import AugmentNoise

        aug = AugmentNoise(sigma=0.05)
        _, ya = aug(self.X, self.y)
        np.testing.assert_array_equal(ya, self.y)

    def test_zero_sigma_identity(self):
        from pycsamt.ai.training.augment import AugmentNoise

        aug = AugmentNoise(sigma=0.0, clip=None)
        Xa, _ = aug(self.X, self.y)
        np.testing.assert_array_equal(Xa, self.X)

    def test_y_none(self):
        from pycsamt.ai.training.augment import AugmentNoise

        aug = AugmentNoise()
        Xa, ya = aug(self.X)
        self.assertIsNone(ya)


class TestAugmentStaticShift(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(1)
        self.X = self.rng.standard_normal((30, 20)).astype(np.float32)

    def test_tuple_range(self):
        from pycsamt.ai.training.augment import (
            AugmentStaticShift,
        )

        aug = AugmentStaticShift(shift_range=(0.3, 3.0))
        Xa, _ = aug(self.X)
        self.assertEqual(Xa.shape, self.X.shape)

    def test_scalar_range_coerced(self):
        """Scalar 0.1 → (10^-0.1, 10^+0.1)."""
        from pycsamt.ai.training.augment import (
            AugmentStaticShift,
        )

        aug = AugmentStaticShift(0.1)
        self.assertAlmostEqual(aug.shift_range[0], 10**-0.1, places=6)
        self.assertAlmostEqual(aug.shift_range[1], 10**0.1, places=6)

    def test_scalar_no_error(self):
        from pycsamt.ai.training.augment import (
            AugmentStaticShift,
        )

        aug = AugmentStaticShift(0.3)
        Xa, _ = aug(self.X)
        self.assertEqual(Xa.shape, self.X.shape)

    def test_per_sample_variation(self):
        from pycsamt.ai.training.augment import (
            AugmentStaticShift,
        )

        aug = AugmentStaticShift(shift_range=(0.01, 100.0), per_sample=True)
        Xa, _ = aug(self.X)
        diffs = np.abs(Xa - self.X).mean(axis=1)
        self.assertGreater(diffs.std(), 0)


class TestAugmentFreqDrop(unittest.TestCase):
    def test_zeros_introduced(self):
        from pycsamt.ai.training.augment import (
            AugmentFreqDrop,
        )

        X = np.ones((20, 40), dtype=np.float32)
        aug = AugmentFreqDrop(drop_rate=0.25)
        Xa, _ = aug(X)
        self.assertGreater((Xa == 0.0).sum(), 0)

    def test_contiguous(self):
        from pycsamt.ai.training.augment import (
            AugmentFreqDrop,
        )

        X = np.ones((10, 50), dtype=np.float32)
        aug = AugmentFreqDrop(drop_rate=0.1, contiguous=True)
        Xa, _ = aug(X)
        self.assertEqual(Xa.shape, X.shape)


class TestAugmentMixup(unittest.TestCase):
    def test_output_in_convex_hull(self):
        from pycsamt.ai.training.augment import AugmentMixup

        X = np.vstack([np.zeros((25, 10)), np.ones((25, 10))]).astype(np.float32)
        y = X.copy()
        aug = AugmentMixup(alpha=1.0)
        Xa, ya = aug(X, y)
        self.assertTrue(np.all(Xa >= 0) and np.all(Xa <= 1))
        self.assertTrue(np.all(ya >= 0) and np.all(ya <= 1))

    def test_single_sample_passthrough(self):
        from pycsamt.ai.training.augment import AugmentMixup

        X = np.ones((1, 10), dtype=np.float32)
        aug = AugmentMixup()
        Xa, _ = aug(X)
        np.testing.assert_array_equal(Xa, X)


class TestCompose(unittest.TestCase):
    def setUp(self):
        from pycsamt.ai.training.augment import (
            AugmentNoise,
            AugmentStaticShift,
            Compose,
        )

        self.aug = Compose(
            [
                AugmentNoise(sigma=0.03),
                AugmentStaticShift(0.2),
            ],
            seed=0,
        )
        self.X = np.random.randn(40, 20).astype(np.float32)
        self.y = np.random.randn(40, 5).astype(np.float32)

    def test_output_shape(self):
        Xa, ya = self.aug(self.X, self.y)
        self.assertEqual(Xa.shape, self.X.shape)
        self.assertEqual(ya.shape, self.y.shape)

    def test_modifies_X(self):
        Xa, _ = self.aug(self.X, self.y)
        self.assertFalse(np.allclose(Xa, self.X))

    def test_y_none(self):
        Xa, ya = self.aug(self.X)
        self.assertIsNone(ya)


class TestRandomApply(unittest.TestCase):
    def test_p0_never_applies(self):
        from pycsamt.ai.training.augment import (
            AugmentNoise,
            RandomApply,
        )

        aug = RandomApply(AugmentNoise(sigma=10.0), p=0.0, seed=0)
        X = np.ones((20, 10), dtype=np.float32)
        for _ in range(10):
            Xa, _ = aug(X)
            np.testing.assert_array_equal(Xa, X)

    def test_p1_always_applies(self):
        from pycsamt.ai.training.augment import (
            AugmentNoise,
            RandomApply,
        )

        aug = RandomApply(AugmentNoise(sigma=5.0), p=1.0, seed=0)
        X = np.ones((20, 10), dtype=np.float32)
        Xa, _ = aug(X)
        self.assertFalse(np.allclose(Xa, X))


# ─────────────────────────────────────────────────────────────────────────────
# Metrics
# ─────────────────────────────────────────────────────────────────────────────


class TestMetrics(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(99)
        self.y_true = self.rng.standard_normal((50, 8)).astype(np.float32)
        self.y_pred = self.y_true + 0.01 * self.rng.standard_normal((50, 8)).astype(
            np.float32
        )

    def test_rmse_near_zero(self):
        from pycsamt.ai.training.metrics import rmse

        val = rmse(self.y_true, self.y_pred)
        self.assertGreater(val, 0)
        self.assertLess(val, 0.1)

    def test_mae_near_zero(self):
        from pycsamt.ai.training.metrics import mae

        val = mae(self.y_true, self.y_pred)
        self.assertGreater(val, 0)
        self.assertLess(val, 0.1)

    def test_r2_near_one(self):
        from pycsamt.ai.training.metrics import r2

        val = r2(self.y_true, self.y_pred)
        self.assertGreater(val, 0.99)

    def test_relative_rmse(self):
        from pycsamt.ai.training.metrics import relative_rmse

        val = relative_rmse(self.y_true, self.y_pred)
        self.assertIsInstance(val, float)

    def test_layer_rmse_shape(self):
        from pycsamt.ai.training.metrics import layer_rmse

        lrmse = layer_rmse(self.y_true, self.y_pred)
        self.assertEqual(lrmse.shape, (8,))

    def test_summarise_keys(self):
        from pycsamt.ai.training.metrics import summarise

        m = summarise(self.y_true, self.y_pred)
        for key in ("rmse", "mae", "r2", "relative_rmse"):
            self.assertIn(key, m)

    def test_depth_rmse(self):
        from pycsamt.ai.training.metrics import depth_rmse

        val = depth_rmse(self.y_true, self.y_pred, n_layers=4)
        self.assertIsInstance(val, float)
        self.assertGreater(val, 0)

    def test_nan_ignored(self):
        from pycsamt.ai.training.metrics import rmse

        yt = self.y_true.copy()
        yt[0, 0] = np.nan
        yp = self.y_pred.copy()
        val = rmse(yt, yp)
        self.assertTrue(np.isfinite(val))

    def test_add_noise_forward(self):
        """add_noise in forward.noise operates on ForwardResponse objects."""
        from pycsamt.forward import LayeredModel, MT1DForward
        from pycsamt.forward.noise import add_noise

        model = LayeredModel(resistivity=[100.0, 10.0, 500.0], thickness=[30.0, 200.0])
        freqs = np.logspace(0, 3, 20)
        resp = MT1DForward(freqs).run(model)
        resp_noisy = add_noise(resp, noise_model="gaussian", level=0.1, seed=0)
        self.assertIsNotNone(resp_noisy.rho_a)
        # Noisy response should differ from clean
        self.assertFalse(np.allclose(resp_noisy.rho_a, resp.rho_a))


if __name__ == "__main__":
    unittest.main()
