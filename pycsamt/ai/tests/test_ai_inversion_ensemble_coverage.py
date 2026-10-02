# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Coverage tests for pycsamt.ai.inversion.ensemble.EnsembleInverter.

Complements the calibration-focused ``TestEnsembleInverterCalibrate`` in
test_inversion.py (which stubs out ``_members`` entirely) and the single
``test_ensemble_predict_shape`` smoke test: this file drives a real
``EnsembleInverter`` through a full ``fit`` on a tiny
:class:`~pycsamt.ai.inversion.inv1d.EMInverter1D` ensemble, then exercises
``predict``/``predict_quantiles``/``score``/``coverage``, the
directory-based ``save``/``load`` round trip, ``plot_uncertainty_profile``,
and the dunder methods.
"""

from __future__ import annotations

import unittest

import numpy as np


def _has_backend():
    from pycsamt.backends import get_backend

    return get_backend() != "none"


@unittest.skipUnless(_has_backend(), "no DL backend available")
class TestEnsembleInverterFitPredict(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pycsamt.forward.batch import generate_dataset

        cls.ds = generate_dataset(
            solver="MT1D",
            n_samples=48,
            freqs=np.logspace(0, 3, 12),
            n_layers=3,
            seed=7,
            verbose=False,
        )

    def _new_ensemble(self, n_estimators=2, **kwargs):
        from pycsamt.ai.inversion import EnsembleInverter
        from pycsamt.ai.inversion.inv1d import EMInverter1D

        base = EMInverter1D(arch="cnn1d", n_layers=3, solver="mt1d")
        return EnsembleInverter(
            base_estimator=base, n_estimators=n_estimators, **kwargs
        )

    def test_seeds_padded_when_fewer_than_n_estimators(self):
        ens = self._new_ensemble(n_estimators=4, seeds=[10, 11])
        self.assertEqual(len(ens.seeds), 4)
        self.assertEqual(ens.seeds[:2], [10, 11])

    def test_repr_unfitted(self):
        ens = self._new_ensemble()
        self.assertIn("unfitted", repr(ens))
        self.assertIn("n_estimators=2", repr(ens))

    def test_fit_returns_self_and_marks_fitted(self):
        ens = self._new_ensemble()
        out = ens.fit(self.ds, epochs=2, batch_size=16, verbose=False)
        self.assertIs(out, ens)
        self.assertTrue(ens._is_fitted)
        self.assertEqual(len(ens), 2)
        self.assertIn("fitted", repr(ens))

    def test_getitem_returns_member(self):
        ens = self._new_ensemble()
        ens.fit(self.ds, epochs=2, batch_size=16, verbose=False)
        member = ens[0]
        self.assertTrue(hasattr(member, "predict"))

    def test_predict_mean_shape(self):
        ens = self._new_ensemble()
        ens.fit(self.ds, epochs=2, batch_size=16, verbose=False)
        y_mean = ens.predict(self.ds.X[:8])
        self.assertEqual(y_mean.shape, (8, self.ds.y.shape[1]))

    def test_predict_quantiles_keys_and_shapes(self):
        ens = self._new_ensemble()
        ens.fit(self.ds, epochs=2, batch_size=16, verbose=False)
        q = ens.predict_quantiles(self.ds.X[:6], q=(0.1, 0.5, 0.9))
        self.assertEqual(set(q.keys()), {0.1, 0.5, 0.9})
        for arr in q.values():
            self.assertEqual(arr.shape, (6, self.ds.y.shape[1]))

    def test_score_returns_finite_float(self):
        ens = self._new_ensemble()
        ens.fit(self.ds, epochs=2, batch_size=16, verbose=False)
        score = ens.score(self.ds.X, self.ds.y, metric="rmse")
        self.assertTrue(np.isfinite(score))

    def test_coverage_within_zero_one(self):
        ens = self._new_ensemble()
        ens.fit(self.ds, epochs=2, batch_size=16, verbose=False)
        cov = ens.coverage(self.ds.X, self.ds.y, n_sigma=1.96)
        self.assertTrue(0.0 <= cov <= 1.0)

    def test_save_load_round_trip(self):
        import tempfile
        from pathlib import Path

        from pycsamt.ai.inversion import EnsembleInverter

        ens = self._new_ensemble()
        ens.fit(self.ds, epochs=2, batch_size=16, verbose=False)
        y_before = ens.predict(self.ds.X[:5])

        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "ens_dir"
            ens.save(path)
            self.assertTrue((path / "member_00.npz").exists())
            self.assertTrue((path / "member_01.npz").exists())
            loaded = EnsembleInverter.load(path)

        self.assertEqual(loaded.n_estimators, 2)
        self.assertTrue(loaded._is_fitted)
        y_after = loaded.predict(self.ds.X[:5])
        np.testing.assert_allclose(y_before, y_after, rtol=1e-4, atol=1e-5)

    def test_plot_uncertainty_profile_returns_figure(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ens = self._new_ensemble()
        ens.fit(self.ds, epochs=2, batch_size=16, verbose=False)
        fig = ens.plot_uncertainty_profile(
            self.ds.X[:3], sample_idx=0, y_true=self.ds.y[0]
        )
        self.assertEqual(fig.__class__.__name__, "Figure")
        plt.close(fig)

    def test_predict_before_fit_raises(self):
        ens = self._new_ensemble()
        with self.assertRaises(RuntimeError):
            ens.predict(self.ds.X[:2])


if __name__ == "__main__":
    unittest.main()
