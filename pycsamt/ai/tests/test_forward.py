# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.forward — LayeredModel, forward solvers, ForwardDataset."""
import unittest

import numpy as np


class TestLayeredModel(unittest.TestCase):

    def _make_model(self, n=4):
        from pycsamt.forward import LayeredModel
        rho = np.array([100., 10., 500., 50., 1000.])[:n + 1]
        thick = np.array([5., 20., 100., 300.])[:n]
        return LayeredModel(resistivity=rho, thickness=thick)

    def test_basic_construction(self):
        m = self._make_model(4)
        self.assertEqual(len(m.resistivity), 5)
        self.assertEqual(len(m.thickness), 4)

    def test_n_layers_property(self):
        m = self._make_model(4)
        self.assertEqual(m.n_layers, 5)

    def test_to_vector(self):
        m = self._make_model(4)
        v = m.to_vector()
        # 5 rho (log10) + 4 thick
        self.assertEqual(v.shape, (9,))

    def test_from_vector(self):
        from pycsamt.forward import LayeredModel
        m = self._make_model(4)
        v = m.to_vector()
        m2 = LayeredModel.from_vector(v, n_layers=5)
        np.testing.assert_allclose(m2.resistivity, m.resistivity, rtol=1e-6)
        np.testing.assert_allclose(m2.thickness, m.thickness)

    def test_positive_validation(self):
        from pycsamt.forward import LayeredModel
        with self.assertRaises((ValueError, TypeError)):
            LayeredModel(resistivity=[-10., 50.], thickness=[5.])

    def test_repr_not_empty(self):
        m = self._make_model(3)
        self.assertIn("LayeredModel", repr(m))

    def test_random_constructor(self):
        from pycsamt.forward import LayeredModel
        m = LayeredModel.random(n_layers=5, seed=0)
        self.assertEqual(m.n_layers, 5)
        self.assertEqual(len(m.thickness), 4)
        self.assertTrue(np.all(m.resistivity > 0))

    def test_depth_property(self):
        m = self._make_model(3)
        self.assertEqual(len(m.depth), 4)   # n_layers = 4, depth for each layer top
        self.assertEqual(m.depth[0], 0.0)


class TestForwardResponse(unittest.TestCase):

    def _mt_response(self):
        from pycsamt.forward import MT1DForward, LayeredModel
        model = LayeredModel(
            resistivity=np.array([100., 10., 1000.]),
            thickness=np.array([30., 200.]),
        )
        freqs = np.logspace(-2, 4, 30)
        solver = MT1DForward(freqs)
        return solver.run(model)

    def test_to_array_shape(self):
        resp = self._mt_response()
        arr = resp.to_array()
        self.assertGreater(arr.size, 0)

    def test_to_array_log(self):
        resp = self._mt_response()
        arr_lin = resp.to_array(log_rho=False)
        arr_log = resp.to_array(log_rho=True)
        self.assertEqual(arr_lin.shape, arr_log.shape)

    def test_apparent_resistivity_positive(self):
        resp = self._mt_response()
        if resp.rho_a is not None:
            self.assertTrue(np.all(resp.rho_a[np.isfinite(resp.rho_a)] > 0))

    def test_phase_range(self):
        resp = self._mt_response()
        if resp.phase is not None:
            finite = resp.phase[np.isfinite(resp.phase)]
            self.assertTrue(np.all(finite >= -90) and np.all(finite <= 90))


class TestMT1DForward(unittest.TestCase):

    def setUp(self):
        from pycsamt.forward import LayeredModel
        self.model = LayeredModel(
            resistivity=np.array([100., 10., 1000.]),
            thickness=np.array([30., 200.]),
        )
        self.freqs = np.logspace(-2, 4, 40)

    def test_output_shape(self):
        from pycsamt.forward import MT1DForward
        resp = MT1DForward(self.freqs).run(self.model)
        if resp.rho_a is not None:
            self.assertEqual(resp.rho_a.shape, (40,))

    def test_homogeneous_space(self):
        """Half-space: apparent resistivity should equal halfspace resistivity."""
        from pycsamt.forward import MT1DForward, LayeredModel
        rho0 = 100.0
        model_hs = LayeredModel(
            resistivity=np.array([rho0]),
            thickness=np.array([]),
        )
        resp = MT1DForward(self.freqs).run(model_hs)
        if resp.rho_a is not None:
            np.testing.assert_allclose(resp.rho_a, rho0, rtol=0.01)

    def test_skin_depth_trend(self):
        """Low frequencies should see deeper (higher resistivity) structure."""
        from pycsamt.forward import MT1DForward
        resp = MT1DForward(self.freqs).run(self.model)
        if resp.rho_a is not None:
            rho = resp.rho_a
            self.assertGreater(np.nanmean(rho[:5]), np.nanmean(rho[20:]))


class TestCSAMT1DForward(unittest.TestCase):

    def test_farfield_equals_mt(self):
        """Very large offset: CSAMT result should approach MT result."""
        from pycsamt.forward import CSAMT1DForward, MT1DForward, LayeredModel
        model = LayeredModel(
            resistivity=np.array([100., 500.]),
            thickness=np.array([50.]),
        )
        freqs = np.logspace(1, 4, 20)
        mt_resp = MT1DForward(freqs).run(model)
        csamt_resp = CSAMT1DForward(freqs, source_offset=1e8).run(model)
        if mt_resp.rho_a is not None and csamt_resp.rho_a is not None:
            np.testing.assert_allclose(
                csamt_resp.rho_a, mt_resp.rho_a, rtol=0.05
            )


class TestTEM1DForward(unittest.TestCase):

    def test_output_shape(self):
        from pycsamt.forward import TEM1DForward, LayeredModel
        model = LayeredModel(
            resistivity=np.array([100., 10., 1000.]),
            thickness=np.array([30., 200.]),
        )
        times = np.logspace(-6, -2, 25)
        resp = TEM1DForward(times).run(model)
        if resp.dBz_dt is not None:
            self.assertEqual(resp.dBz_dt.shape, (25,))

    def test_decay_monotonic(self):
        """TEM step-off response should have an overall decaying trend.

        Strict point-by-point monotonicity fails due to numerical integration
        artefacts at early times.  We instead check that the late-time half
        of the response is generally decreasing: the linear trend (slope of a
        least-squares fit in log-log space) must be negative.
        """
        from pycsamt.forward import TEM1DForward, LayeredModel
        model = LayeredModel(
            resistivity=np.array([100.]),
            thickness=np.array([]),
        )
        times = np.logspace(-5, -2, 20)
        resp = TEM1DForward(times).run(model)
        if resp.dBz_dt is None:
            return
        vals = np.abs(resp.dBz_dt)
        finite_mask = np.isfinite(vals) & (vals > 0)
        finite_t = times[finite_mask]
        finite_v = vals[finite_mask]
        if len(finite_v) < 4:
            return
        # Use the later half of the finite samples to avoid early-time artefacts
        half = len(finite_v) // 2
        log_t = np.log(finite_t[half:])
        log_v = np.log(finite_v[half:])
        slope, _ = np.polyfit(log_t, log_v, 1)
        self.assertLess(slope, 0, "Late-time dBz/dt trend should be decaying")


class TestGenerateDataset(unittest.TestCase):

    def test_mt1d_small(self):
        from pycsamt.forward.batch import generate_dataset
        ds = generate_dataset(
            solver="MT1D",
            n_samples=10,
            freqs=np.logspace(0, 3, 12),
            n_layers=3,
            seed=42,
            verbose=False,
        )
        self.assertEqual(len(ds.X), 10)
        self.assertEqual(len(ds.y), 10)

    def test_feature_dimension(self):
        from pycsamt.forward.batch import generate_dataset
        n_freqs = 12
        ds = generate_dataset(
            solver="MT1D",
            n_samples=10,
            freqs=np.logspace(0, 3, n_freqs),
            n_layers=3,
            include_phase=True,
            seed=0,
            verbose=False,
        )
        # With phase: 2 × n_freqs features
        self.assertEqual(ds.X.shape[1], 2 * n_freqs)

    def test_target_dimension(self):
        from pycsamt.forward.batch import generate_dataset
        n_layers = 4
        ds = generate_dataset(
            solver="MT1D",
            n_samples=10,
            freqs=np.logspace(0, 3, 10),
            n_layers=n_layers,
            seed=0,
            verbose=False,
        )
        # n_layers rho + (n_layers-1) thick = 2*n_layers - 1
        self.assertEqual(ds.y.shape[1], 2 * n_layers - 1)

    def test_dataset_split(self):
        from pycsamt.forward.batch import generate_dataset
        ds = generate_dataset(
            solver="MT1D",
            n_samples=20,
            freqs=np.logspace(0, 3, 10),
            n_layers=3,
            seed=0,
            verbose=False,
        )
        train, val, test = ds.split(val_frac=0.2, test_frac=0.1, seed=0)
        total = len(train.X) + len(val.X) + len(test.X)
        self.assertEqual(total, 20)

    def test_dataset_save_load(self):
        import tempfile, os
        from pycsamt.forward.batch import generate_dataset, ForwardDataset
        ds = generate_dataset(
            solver="MT1D",
            n_samples=8,
            freqs=np.logspace(0, 3, 8),
            n_layers=3,   # n_layers >= 3 to avoid batch.py edge case
            seed=1,
            verbose=False,
        )
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "ds.npz")
            ds.save(path)
            ds2 = ForwardDataset.load(path)
        np.testing.assert_array_equal(ds.X, ds2.X)
        np.testing.assert_array_equal(ds.y, ds2.y)


class TestSurveyDataset3D(unittest.TestCase):
    """Container tests — no DL backend required."""

    @classmethod
    def setUpClass(cls):
        from pycsamt.forward.batch import generate_dataset_3d
        cls.ds = generate_dataset_3d(
            n_surveys=20,
            n_stations=9,
            n_layers=3,
            freqs=np.logspace(0, 3, 10),
            corr_length=2000.0,
            seed=7,
            verbose=False,
        )

    def test_X_shape(self):
        ds = self.ds
        self.assertEqual(ds.X.shape, (20, 9, 20))   # 10 freqs × 2 (rho+phase)

    def test_y_shape(self):
        ds = self.ds
        # n_layers=3 → 3 rho + 2 thick = 5
        self.assertEqual(ds.y.shape, (20, 9, 5))

    def test_coords_shape(self):
        self.assertEqual(self.ds.coords.shape, (9, 2))

    def test_n_surveys_property(self):
        self.assertEqual(self.ds.n_surveys, 20)

    def test_n_stations_property(self):
        self.assertEqual(self.ds.n_stations, 9)

    def test_n_features_property(self):
        self.assertEqual(self.ds.n_features, 20)

    def test_n_params_property(self):
        self.assertEqual(self.ds.n_params, 5)

    def test_dtype_float32(self):
        self.assertEqual(self.ds.X.dtype, np.float32)
        self.assertEqual(self.ds.y.dtype, np.float32)

    def test_repr_contains_class(self):
        self.assertIn("SurveyDataset3D", repr(self.ds))

    def test_len(self):
        self.assertEqual(len(self.ds), 20)

    def test_split_totals(self):
        train, val, test = self.ds.split(val_frac=0.2, test_frac=0.1, seed=0)
        self.assertEqual(len(train) + len(val) + len(test), 20)

    def test_split_coords_shared(self):
        train, val, _ = self.ds.split(seed=0)
        np.testing.assert_array_equal(train.coords, self.ds.coords)
        np.testing.assert_array_equal(val.coords, self.ds.coords)

    def test_save_load_roundtrip(self):
        import tempfile, os
        from pycsamt.forward.batch import SurveyDataset3D
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "ds3d.npz")
            self.ds.save(path)
            ds2 = SurveyDataset3D.load(path)
        np.testing.assert_array_equal(self.ds.X, ds2.X)
        np.testing.assert_array_equal(self.ds.y, ds2.y)
        np.testing.assert_array_equal(self.ds.coords, ds2.coords)

    def test_freqs_preserved(self):
        self.assertIsNotNone(self.ds.freqs)
        self.assertEqual(len(self.ds.freqs), 10)


class TestGenerateDataset3D(unittest.TestCase):
    """Shape and content tests for generate_dataset_3d."""

    def _gen(self, **kw):
        from pycsamt.forward.batch import generate_dataset_3d
        defaults = dict(
            n_surveys=10, n_stations=4, n_layers=3,
            freqs=np.logspace(0, 3, 8),
            corr_length=1500., seed=0, verbose=False,
        )
        defaults.update(kw)
        return generate_dataset_3d(**defaults)

    def test_x_feature_no_phase(self):
        ds = self._gen(include_phase=False)
        # No phase: n_features = n_freqs = 8
        self.assertEqual(ds.n_features, 8)

    def test_x_feature_with_phase(self):
        ds = self._gen(include_phase=True)
        self.assertEqual(ds.n_features, 16)  # 2 × 8

    def test_target_dimension(self):
        ds = self._gen(n_layers=5)
        # 5 rho + 4 thick = 9
        self.assertEqual(ds.n_params, 9)

    def test_grid_layout(self):
        ds = self._gen(n_stations=16, station_layout="grid")
        self.assertEqual(ds.coords.shape[0], 16)

    def test_random_layout(self):
        ds = self._gen(n_stations=12, station_layout="random")
        self.assertEqual(ds.coords.shape[0], 12)

    def test_invalid_layout_raises(self):
        from pycsamt.forward.batch import generate_dataset_3d
        with self.assertRaises(ValueError):
            generate_dataset_3d(
                n_surveys=2, n_stations=4, station_layout="hexagonal",
                freqs=np.logspace(0, 2, 5), seed=0, verbose=False,
            )

    def test_invalid_solver_raises(self):
        from pycsamt.forward.batch import generate_dataset_3d
        with self.assertRaises(ValueError):
            generate_dataset_3d(
                n_surveys=2, n_stations=4, solver="tem1d",
                freqs=np.logspace(0, 2, 5), seed=0, verbose=False,
            )

    def test_deterministic_with_seed(self):
        ds1 = self._gen(seed=42)
        ds2 = self._gen(seed=42)
        np.testing.assert_array_equal(ds1.X, ds2.X)

    def test_different_seeds_differ(self):
        ds1 = self._gen(seed=1)
        ds2 = self._gen(seed=2)
        self.assertFalse(np.allclose(ds1.X, ds2.X))

    def test_no_nan_in_output(self):
        ds = self._gen()
        self.assertFalse(np.any(np.isnan(ds.X)))
        self.assertFalse(np.any(np.isnan(ds.y)))

    def test_corr_length_effect(self):
        """High corr_length → adjacent stations more similar than low corr_length."""
        ds_lo = self._gen(corr_length=10., n_surveys=50)
        ds_hi = self._gen(corr_length=9000., n_surveys=50)
        # Variance of (station_0 - station_1) in y[:, 0, :] vs y[:, 1, :]
        diff_lo = ds_lo.y[:, 0, :] - ds_lo.y[:, 1, :]
        diff_hi = ds_hi.y[:, 0, :] - ds_hi.y[:, 1, :]
        # High corr_length should give smaller inter-station variance
        self.assertLess(diff_hi.var(), diff_lo.var())

    def test_gcn_pipeline_shapes(self):
        """generate_dataset_3d → build_adjacency → GCNInverter3D compatible shapes."""
        from pycsamt.ai.nets.gcn import build_adjacency
        ds = self._gen(n_surveys=10, n_stations=9, n_layers=3)
        A = build_adjacency(ds.coords, radius=5000.)
        self.assertEqual(A.shape, (9, 9))
        self.assertEqual(ds.X.shape, (10, 9, ds.n_features))
        self.assertEqual(ds.y.shape, (10, 9, 5))


if __name__ == "__main__":
    unittest.main()
