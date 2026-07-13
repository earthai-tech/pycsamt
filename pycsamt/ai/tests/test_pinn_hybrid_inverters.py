# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""Tests for PINNInverter[1-3]D and HybridInverter[1-3]D.

Backend matrix
--------------
Run with the ``pycsamt-v2`` conda environment for the
PyTorch test classes (``TestPINNTorch``,
``TestHybridTorch``) and with the ``base-attentive``
environment for the TensorFlow classes
(``TestPINNTF``, ``TestHybridTF``).
Skip decorators ensure only the tests whose backend
is installed actually execute.
"""

import unittest

import numpy as np

# ── global constants ─────────────────────────────────

_NF = 8  # n_freqs (small for speed)
_NS = 3  # n_stations
_NL = 3  # n_layers
_FREQS = np.logspace(3, 0, _NF)  # 1000 … 1 Hz


# ── backend probe helpers ─────────────────────────────


def _has_torch():
    try:
        import torch  # noqa: F401

        return True
    except ImportError:
        return False


def _has_tf():
    try:
        import tensorflow  # noqa: F401

        return True
    except ImportError:
        return False


def _set_backend(name):
    from pycsamt.backends import set_backend

    set_backend(name)


# ── SiteObs factories ─────────────────────────────────


def _obs1d(name="S01"):
    from pycsamt.ai.inversion._sites_bridge import (
        SiteObs1D,
    )

    rho = np.full(_NF, 100.0)
    ph = np.full(_NF, 45.0)
    return SiteObs1D(
        name=name,
        freq=_FREQS.copy(),
        rho_obs=rho,
        phase_obs=ph,
    )


def _obs2d(name="S01"):
    from pycsamt.ai.inversion._sites_bridge import (
        SiteObs2D,
    )

    rho = np.full(_NF, 100.0)
    ph = np.full(_NF, 45.0)
    return SiteObs2D(
        name=name,
        freq=_FREQS.copy(),
        rho_te=rho,
        phase_te=ph,
        rho_tm=rho.copy(),
        phase_tm=ph.copy(),
    )


def _obs_list_1d(n=_NS):
    return [_obs1d(f"S{i + 1:02d}") for i in range(n)]


def _obs_list_2d(n=_NS):
    return [_obs2d(f"S{i + 1:02d}") for i in range(n)]


def _freqs_grid_for(obs_list):
    from pycsamt.ai.inversion._sites_bridge import (
        _make_common_grid,
    )

    return _make_common_grid(
        [o.freq for o in obs_list],
        n_freqs=_NF,
    )


def _adj(n=_NS):
    A = np.ones((n, n)) - np.eye(n)
    return A.astype(np.float64)


def _coords(n=_NS):
    return np.column_stack([np.arange(n) * 500.0, np.zeros(n)])


# ── mock AI inverters ─────────────────────────────────


class _MockAI1D:
    _is_fitted = True
    n_layers = _NL
    _n_features = _NF * 2

    def predict(self, X, as_log_rho=False):
        n = X.shape[0]
        return np.full((n, 2 * _NL - 1), 2.0)


class _MockAI2D:
    _is_fitted = True
    n_depth = _NL
    n_components = 2
    n_freqs = _NF
    n_stations = _NS

    def predict(self, panel, as_log_rho=False):
        n_sta = panel.shape[3]
        return np.full((1, self.n_depth, n_sta), 2.0)


class _MockAI3D:
    _is_fitted = True
    n_layers = _NL
    n_features = _NF * 2

    def predict(
        self,
        X,
        adjacency=None,
        as_log_rho=False,
    ):
        n = X.shape[0]
        return np.full((n, 2 * _NL - 1), 2.0)


# ── inverter factories (bypass __init__) ─────────────


def _make_pinn1d(**kw):
    from pycsamt.ai.inversion import PINNInverter1D

    inv = object.__new__(PINNInverter1D)
    inv.n_layers = kw.get("n_layers", _NL)
    inv.depth_max = kw.get("depth_max", 500.0)
    inv.device = kw.get("device", None)
    inv._is_fitted = False
    inv._history = []
    inv.solver = "mt1d"
    inv.smoothness_weight = 0.01
    inv.lr = 0.01
    inv.comp = "xy"
    inv.recursive = True
    inv.on_dup = "replace"
    inv.verbose = 0
    inv._obs = _obs_list_1d(kw.get("n_stations", _NS))
    inv._results = []
    return inv


def _make_pinn2d(**kw):
    from pycsamt.ai.inversion import PINNInverter2D

    n_s = kw.get("n_stations", _NS)
    obs = _obs_list_2d(n_s)
    inv = object.__new__(PINNInverter2D)
    inv.n_layers = kw.get("n_layers", _NL)
    inv.depth_max = kw.get("depth_max", 500.0)
    inv.device = kw.get("device", None)
    inv._is_fitted = False
    inv._history = []
    inv.n_freqs = _NF
    inv.mode = "te"
    inv.smoothness_weight = 0.01
    inv.lateral_weight = 0.005
    inv.epochs = kw.get("epochs", 5)
    inv.lr = 0.01
    inv.comp_te = "xy"
    inv.comp_tm = "yx"
    inv.verbose = 0
    inv._obs = obs
    inv._freqs_grid = _freqs_grid_for(obs)
    inv._result = None
    return inv


def _make_pinn3d(**kw):
    from pycsamt.ai.inversion import PINNInverter3D

    n_s = kw.get("n_stations", _NS)
    obs = _obs_list_2d(n_s)
    inv = object.__new__(PINNInverter3D)
    inv.n_layers = kw.get("n_layers", _NL)
    inv.depth_max = kw.get("depth_max", 500.0)
    inv.device = kw.get("device", None)
    inv._is_fitted = False
    inv._history = []
    inv.n_freqs = _NF
    inv.mode = "te"
    inv.smoothness_weight = 0.01
    inv.graph_weight = 0.003
    inv.radius = 5000.0
    inv.epochs = kw.get("epochs", 5)
    inv.lr = 0.01
    inv.comp_te = "xy"
    inv.comp_tm = "yx"
    inv.verbose = 0
    inv._obs = obs
    inv._freqs_grid = _freqs_grid_for(obs)
    inv._adjacency = _adj(n_s)
    inv._coords = _coords(n_s)
    inv._result = None
    return inv


def _make_hybrid1d(**kw):
    from pycsamt.ai.inversion import HybridInverter1D

    n_s = kw.get("n_stations", _NS)
    inv = object.__new__(HybridInverter1D)
    inv.n_layers = _NL
    inv.depth_max = 500.0
    inv.device = kw.get("device", None)
    inv._is_fitted = False
    inv._history = []
    inv._ai_inv = _MockAI1D()
    inv._stage1_fitted = False
    inv.solver = "mt1d"
    inv.max_iter = kw.get("max_iter", 5)
    inv.smoothness_weight = 0.005
    inv.lr = 5e-3
    inv.comp = "xy"
    inv.n_freqs = _NF
    inv.recursive = True
    inv.on_dup = "replace"
    inv.verbose = 0
    inv._obs = _obs_list_1d(n_s)
    inv._stage1 = []
    inv._stage2 = []
    return inv


def _make_hybrid2d(**kw):
    from pycsamt.ai.inversion import HybridInverter2D

    n_s = kw.get("n_stations", _NS)
    obs = _obs_list_2d(n_s)
    inv = object.__new__(HybridInverter2D)
    inv.n_layers = _NL
    inv.depth_max = kw.get("depth_max", 500.0)
    inv.device = kw.get("device", None)
    inv._is_fitted = False
    inv._history = []
    inv._ai_inv = _MockAI2D()
    inv._stage1_fitted = False
    inv.n_freqs = _NF
    inv.mode = "te"
    inv.smoothness_weight = 0.005
    inv.lateral_weight = 0.003
    inv.epochs = kw.get("epochs", 5)
    inv.lr = 5e-3
    inv.comp_te = "xy"
    inv.comp_tm = "yx"
    inv.verbose = 0
    inv._obs = obs
    inv._freqs_grid = _freqs_grid_for(obs)
    inv._stage1_log_rho = None
    inv._result = None
    return inv


def _make_hybrid3d(**kw):
    from pycsamt.ai.inversion import HybridInverter3D

    n_s = kw.get("n_stations", _NS)
    obs = _obs_list_2d(n_s)
    inv = object.__new__(HybridInverter3D)
    inv.n_layers = _NL
    inv.depth_max = kw.get("depth_max", 500.0)
    inv.device = kw.get("device", None)
    inv._is_fitted = False
    inv._history = []
    inv._ai_inv = _MockAI3D()
    inv._stage1_fitted = False
    inv.n_freqs = _NF
    inv.mode = "te"
    inv.smoothness_weight = 0.005
    inv.graph_weight = 0.003
    inv.radius = 5000.0
    inv.epochs = kw.get("epochs", 5)
    inv.lr = 5e-3
    inv.comp_te = "xy"
    inv.comp_tm = "yx"
    inv.verbose = 0
    inv._obs = obs
    inv._freqs_grid = _freqs_grid_for(obs)
    inv._adjacency = _adj(n_s)
    inv._coords = _coords(n_s)
    inv._stage1_params = None
    inv._result = None
    return inv


# ════════════════════════════════════════════════════
# 1.  Interface tests (no DL backend required)
# ════════════════════════════════════════════════════


class TestPINNHybridInterface(unittest.TestCase):
    r"""Contract tests; no DL backend required."""

    # ── importability ────────────────────────────────

    def test_pinn_classes_importable(self):
        from pycsamt.ai.inversion import (
            PINNInverter1D,
            PINNInverter2D,
            PINNInverter3D,
        )

        for cls in (
            PINNInverter1D,
            PINNInverter2D,
            PINNInverter3D,
        ):
            self.assertIsNotNone(cls)

    def test_hybrid_classes_importable(self):
        from pycsamt.ai.inversion import (
            HybridInverter1D,
            HybridInverter2D,
            HybridInverter3D,
        )

        for cls in (
            HybridInverter1D,
            HybridInverter2D,
            HybridInverter3D,
        ):
            self.assertIsNotNone(cls)

    # ── class hierarchy ──────────────────────────────

    def test_pinn_inherit_base_pinn(self):
        from pycsamt.ai._base import BasePINNInverter
        from pycsamt.ai.inversion import (
            PINNInverter1D,
            PINNInverter2D,
            PINNInverter3D,
        )

        for cls in (
            PINNInverter1D,
            PINNInverter2D,
            PINNInverter3D,
        ):
            self.assertTrue(
                issubclass(cls, BasePINNInverter),
                msg=cls.__name__,
            )

    def test_hybrid_inherit_base_hybrid(self):
        from pycsamt.ai._base import (
            BaseHybridInverter,
        )
        from pycsamt.ai.inversion import (
            HybridInverter1D,
            HybridInverter2D,
            HybridInverter3D,
        )

        for cls in (
            HybridInverter1D,
            HybridInverter2D,
            HybridInverter3D,
        ):
            self.assertTrue(
                issubclass(cls, BaseHybridInverter),
                msg=cls.__name__,
            )

    def test_all_inherit_pycsamt_object(self):
        from pycsamt.ai.inversion import (
            HybridInverter1D,
            HybridInverter2D,
            HybridInverter3D,
            PINNInverter1D,
            PINNInverter2D,
            PINNInverter3D,
        )
        from pycsamt.api.property import PyCSAMTObject

        for cls in (
            PINNInverter1D,
            PINNInverter2D,
            PINNInverter3D,
            HybridInverter1D,
            HybridInverter2D,
            HybridInverter3D,
        ):
            self.assertTrue(
                issubclass(cls, PyCSAMTObject),
                msg=cls.__name__,
            )

    # ── repr ─────────────────────────────────────────

    def test_repr_pinn1d(self):
        inv = _make_pinn1d()
        self.assertIn("PINNInverter1D", repr(inv))

    def test_repr_pinn2d(self):
        inv = _make_pinn2d()
        self.assertIn("PINNInverter2D", repr(inv))

    def test_repr_pinn3d(self):
        inv = _make_pinn3d()
        self.assertIn("PINNInverter3D", repr(inv))

    def test_repr_hybrid1d(self):
        inv = _make_hybrid1d()
        self.assertIn("HybridInverter1D", repr(inv))

    def test_repr_hybrid2d(self):
        inv = _make_hybrid2d()
        self.assertIn("HybridInverter2D", repr(inv))

    def test_repr_hybrid3d(self):
        inv = _make_hybrid3d()
        self.assertIn("HybridInverter3D", repr(inv))

    # ── unfitted guard ───────────────────────────────

    def test_unfitted_pinn1d_raises(self):
        inv = _make_pinn1d()
        with self.assertRaises(RuntimeError):
            inv.convergence_curve()

    def test_unfitted_pinn2d_raises(self):
        inv = _make_pinn2d()
        with self.assertRaises(RuntimeError):
            inv.convergence_curve()

    def test_unfitted_pinn3d_raises(self):
        inv = _make_pinn3d()
        with self.assertRaises(RuntimeError):
            inv.convergence_curve()

    def test_unfitted_hybrid1d_raises(self):
        inv = _make_hybrid1d()
        with self.assertRaises(RuntimeError):
            inv.predict()

    def test_unfitted_hybrid2d_raises(self):
        inv = _make_hybrid2d()
        with self.assertRaises(RuntimeError):
            inv.resistivity_section()

    def test_unfitted_hybrid3d_raises(self):
        inv = _make_hybrid3d()
        with self.assertRaises(RuntimeError):
            inv.resistivity_volume()

    # ── base methods present ─────────────────────────

    def test_base_methods_on_pinn1d(self):
        inv = _make_pinn1d()
        for meth in (
            "_resolve_device",
            "_require_backend",
            "_check_fitted",
            "convergence_curve",
        ):
            self.assertTrue(
                callable(getattr(inv, meth)),
                msg=meth,
            )

    def test_base_methods_on_hybrid3d(self):
        inv = _make_hybrid3d()
        for meth in (
            "_resolve_device",
            "_require_backend",
            "_check_fitted",
        ):
            self.assertTrue(
                callable(getattr(inv, meth)),
                msg=meth,
            )

    def test_stations_property_pinn2d(self):
        inv = _make_pinn2d()
        names = inv.stations
        self.assertEqual(len(names), _NS)

    def test_n_sites_property_pinn3d(self):
        inv = _make_pinn3d()
        self.assertEqual(inv.n_sites, _NS)


# ════════════════════════════════════════════════════
# 2.  PyTorch backend tests
#     Run with:  conda run -n pycsamt-v2 python -m
#                pytest test_pinn_hybrid_inverters.py
# ════════════════════════════════════════════════════

_SKIP_TORCH = unittest.skipUnless(_has_torch(), "torch not available")


@_SKIP_TORCH
class TestPINNTorch(unittest.TestCase):
    r"""Full fit/output tests against the Torch backend."""

    @classmethod
    def setUpClass(cls):
        _set_backend("torch")

    @classmethod
    def tearDownClass(cls):
        _set_backend("auto")

    # ── PINNInverter1D ───────────────────────────────

    def test_1d_fit_sets_is_fitted(self):
        inv = _make_pinn1d()
        inv.fit(epochs=5, verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_1d_history_row_count(self):
        inv = _make_pinn1d()
        inv.fit(epochs=7, verbose=False)
        df = inv.loss_curves()
        self.assertEqual(len(df), 7 * _NS)

    def test_1d_predict_models_length(self):
        inv = _make_pinn1d()
        inv.fit(epochs=3, verbose=False)
        models = inv.predict()
        self.assertEqual(len(models), _NS)

    def test_1d_residuals_columns(self):
        inv = _make_pinn1d()
        inv.fit(epochs=3, verbose=False)
        df = inv.residuals()
        expected = {
            "station",
            "freq",
            "rho_obs",
            "rho_pred",
            "phase_obs",
            "phase_pred",
        }
        self.assertEqual(set(df.columns), expected)

    def test_1d_convergence_curve_base(self):
        inv = _make_pinn1d()
        inv.fit(epochs=4, verbose=False)
        hist = inv.convergence_curve()
        self.assertIsInstance(hist, list)
        self.assertEqual(len(hist), 0)

    # ── PINNInverter2D ───────────────────────────────

    def test_2d_fit_sets_is_fitted(self):
        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_2d_section_shape(self):
        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        sec = inv.resistivity_section()
        self.assertEqual(sec.shape, (_NL, _NS))

    def test_2d_section_linear(self):
        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        rho = inv.resistivity_section(as_log10=False)
        self.assertTrue(np.all(rho > 0))

    def test_2d_thickness_shape(self):
        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        th = inv.thickness_section()
        exp = (_NL - 1, _NS)
        self.assertEqual(th.shape, exp)

    def test_2d_convergence_curve_df(self):
        import pandas as pd

        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        df = inv.convergence_curve()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn("epoch", df.columns)
        self.assertIn("loss", df.columns)
        self.assertEqual(len(df), 5)

    def test_2d_history_stored(self):
        inv = _make_pinn2d(epochs=6)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._history), 6)

    def test_2d_finite_section(self):
        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        lr = inv.resistivity_section(as_log10=True)
        self.assertTrue(np.all(np.isfinite(lr)))

    # ── PINNInverter3D ───────────────────────────────

    def test_3d_fit_sets_is_fitted(self):
        inv = _make_pinn3d(epochs=5)
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_3d_volume_shape(self):
        inv = _make_pinn3d(epochs=5)
        inv.fit(verbose=False)
        vol = inv.resistivity_volume()
        self.assertEqual(vol.shape, (_NL, _NS))

    def test_3d_volume_linear(self):
        inv = _make_pinn3d(epochs=5)
        inv.fit(verbose=False)
        rho = inv.resistivity_volume(as_log10=False)
        self.assertTrue(np.all(rho > 0))

    def test_3d_thickness_shape(self):
        inv = _make_pinn3d(epochs=5)
        inv.fit(verbose=False)
        th = inv.thickness_volume()
        exp = (_NL - 1, _NS)
        self.assertEqual(th.shape, exp)

    def test_3d_convergence_curve_df(self):
        import pandas as pd

        inv = _make_pinn3d(epochs=5)
        inv.fit(verbose=False)
        df = inv.convergence_curve()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 5)

    def test_3d_history_stored(self):
        inv = _make_pinn3d(epochs=4)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._history), 4)

    def test_3d_adjacency_returns_copy(self):
        inv = _make_pinn3d()
        A = inv.adjacency()
        self.assertEqual(A.shape, (_NS, _NS))
        A[0, 0] = -999.0
        self.assertNotEqual(inv.adjacency()[0, 0], -999.0)

    def test_3d_station_coords_shape(self):
        inv = _make_pinn3d()
        c = inv.station_coords()
        self.assertEqual(c.shape, (_NS, 2))


# ════════════════════════════════════════════════════
# 3.  Hybrid inverter tests — Torch backend
# ════════════════════════════════════════════════════


@_SKIP_TORCH
class TestHybridTorch(unittest.TestCase):
    r"""Hybrid two-stage fit/output tests (Torch)."""

    @classmethod
    def setUpClass(cls):
        _set_backend("torch")

    @classmethod
    def tearDownClass(cls):
        _set_backend("auto")

    # ── HybridInverter1D ─────────────────────────────

    def test_h1d_fit_sets_is_fitted(self):
        inv = _make_hybrid1d(max_iter=5)
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_h1d_stage2_populated(self):
        inv = _make_hybrid1d(max_iter=5)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._stage2), _NS)

    def test_h1d_stage1_populated(self):
        inv = _make_hybrid1d(max_iter=5)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._stage1), _NS)

    def test_h1d_predict_length(self):
        inv = _make_hybrid1d(max_iter=5)
        inv.fit(verbose=False)
        models = inv.predict()
        self.assertEqual(len(models), _NS)

    def test_h1d_stage1_models_length(self):
        inv = _make_hybrid1d(max_iter=5)
        inv.fit(verbose=False)
        m1 = inv.stage1_models()
        self.assertEqual(len(m1), _NS)

    def test_h1d_convergence_curves_df(self):
        import pandas as pd

        inv = _make_hybrid1d(max_iter=5)
        inv.fit(verbose=False)
        df = inv.convergence_curves()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertIn("station", df.columns)
        self.assertIn("loss", df.columns)
        self.assertEqual(len(df), 5 * _NS)

    def test_h1d_stations_property(self):
        inv = _make_hybrid1d(max_iter=3)
        inv.fit(verbose=False)
        self.assertEqual(len(inv.stations), _NS)

    # ── HybridInverter2D ─────────────────────────────

    def test_h2d_fit_sets_is_fitted(self):
        inv = _make_hybrid2d(epochs=5)
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_h2d_section_shape(self):
        inv = _make_hybrid2d(epochs=5)
        inv.fit(verbose=False)
        sec = inv.resistivity_section()
        self.assertEqual(sec.shape, (_NL, _NS))

    def test_h2d_thickness_shape(self):
        inv = _make_hybrid2d(epochs=5)
        inv.fit(verbose=False)
        th = inv.thickness_section()
        exp = (_NL - 1, _NS)
        self.assertEqual(th.shape, exp)

    def test_h2d_stage1_section_shape(self):
        inv = _make_hybrid2d(epochs=5)
        inv.fit(verbose=False)
        s1 = inv.stage1_section()
        self.assertEqual(s1.shape, (_NL, _NS))

    def test_h2d_convergence_curve_df(self):
        import pandas as pd

        inv = _make_hybrid2d(epochs=5)
        inv.fit(verbose=False)
        df = inv.convergence_curve()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 5)

    def test_h2d_history_stored(self):
        inv = _make_hybrid2d(epochs=6)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._history), 6)

    # ── HybridInverter3D ─────────────────────────────

    def test_h3d_fit_sets_is_fitted(self):
        inv = _make_hybrid3d(epochs=5)
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_h3d_volume_shape(self):
        inv = _make_hybrid3d(epochs=5)
        inv.fit(verbose=False)
        vol = inv.resistivity_volume()
        self.assertEqual(vol.shape, (_NL, _NS))

    def test_h3d_thickness_shape(self):
        inv = _make_hybrid3d(epochs=5)
        inv.fit(verbose=False)
        th = inv.thickness_volume()
        exp = (_NL - 1, _NS)
        self.assertEqual(th.shape, exp)

    def test_h3d_stage1_volume_shape(self):
        inv = _make_hybrid3d(epochs=5)
        inv.fit(verbose=False)
        s1 = inv.stage1_volume()
        self.assertEqual(s1.shape, (_NL, _NS))

    def test_h3d_convergence_curve_df(self):
        import pandas as pd

        inv = _make_hybrid3d(epochs=5)
        inv.fit(verbose=False)
        df = inv.convergence_curve()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 5)

    def test_h3d_history_stored(self):
        inv = _make_hybrid3d(epochs=4)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._history), 4)

    def test_h3d_station_coords_shape(self):
        inv = _make_hybrid3d()
        c = inv.station_coords()
        self.assertEqual(c.shape, (_NS, 2))


# ════════════════════════════════════════════════════
# 4.  TensorFlow backend tests
#     Run with:  conda run -n base-attentive python
#                -m pytest test_pinn_hybrid_inverters.py
# ════════════════════════════════════════════════════

_SKIP_TF = unittest.skipUnless(_has_tf(), "tensorflow not available")


@_SKIP_TF
class TestPINNTF(unittest.TestCase):
    r"""Fit/output tests against the TensorFlow backend."""

    @classmethod
    def setUpClass(cls):
        _set_backend("tensorflow")

    @classmethod
    def tearDownClass(cls):
        _set_backend("auto")

    # ── PINNInverter1D ───────────────────────────────

    def test_1d_fit_sets_is_fitted(self):
        inv = _make_pinn1d()
        inv.fit(epochs=5, verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_1d_history_row_count(self):
        inv = _make_pinn1d()
        inv.fit(epochs=6, verbose=False)
        df = inv.loss_curves()
        self.assertEqual(len(df), 6 * _NS)

    def test_1d_residuals_returns_df(self):
        import pandas as pd

        inv = _make_pinn1d()
        inv.fit(epochs=3, verbose=False)
        df = inv.residuals()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertGreater(len(df), 0)

    # ── PINNInverter2D ───────────────────────────────

    def test_2d_fit_sets_is_fitted(self):
        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_2d_section_shape(self):
        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        sec = inv.resistivity_section()
        self.assertEqual(sec.shape, (_NL, _NS))

    def test_2d_section_finite(self):
        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        lr = inv.resistivity_section()
        self.assertTrue(np.all(np.isfinite(lr)))

    def test_2d_thickness_shape(self):
        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        th = inv.thickness_section()
        exp = (_NL - 1, _NS)
        self.assertEqual(th.shape, exp)

    def test_2d_history_length(self):
        inv = _make_pinn2d(epochs=5)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._history), 5)

    # ── PINNInverter3D ───────────────────────────────

    def test_3d_fit_sets_is_fitted(self):
        inv = _make_pinn3d(epochs=5)
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_3d_volume_shape(self):
        inv = _make_pinn3d(epochs=5)
        inv.fit(verbose=False)
        vol = inv.resistivity_volume()
        self.assertEqual(vol.shape, (_NL, _NS))

    def test_3d_thickness_shape(self):
        inv = _make_pinn3d(epochs=5)
        inv.fit(verbose=False)
        th = inv.thickness_volume()
        exp = (_NL - 1, _NS)
        self.assertEqual(th.shape, exp)

    def test_3d_history_length(self):
        inv = _make_pinn3d(epochs=4)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._history), 4)


# ════════════════════════════════════════════════════
# 5.  Hybrid inverter tests — TensorFlow backend
# ════════════════════════════════════════════════════


@_SKIP_TF
class TestHybridTF(unittest.TestCase):
    r"""Hybrid two-stage fit/output tests (TF)."""

    @classmethod
    def setUpClass(cls):
        _set_backend("tensorflow")

    @classmethod
    def tearDownClass(cls):
        _set_backend("auto")

    # ── HybridInverter1D ─────────────────────────────

    def test_h1d_fit_sets_is_fitted(self):
        inv = _make_hybrid1d(max_iter=5)
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_h1d_stage2_populated(self):
        inv = _make_hybrid1d(max_iter=5)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._stage2), _NS)

    def test_h1d_history_row_count(self):
        import pandas as pd

        inv = _make_hybrid1d(max_iter=5)
        inv.fit(verbose=False)
        df = inv.convergence_curves()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 5 * _NS)

    # ── HybridInverter2D ─────────────────────────────

    def test_h2d_fit_sets_is_fitted(self):
        inv = _make_hybrid2d(epochs=5)
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_h2d_section_shape(self):
        inv = _make_hybrid2d(epochs=5)
        inv.fit(verbose=False)
        sec = inv.resistivity_section()
        self.assertEqual(sec.shape, (_NL, _NS))

    def test_h2d_stage1_section_shape(self):
        inv = _make_hybrid2d(epochs=5)
        inv.fit(verbose=False)
        s1 = inv.stage1_section()
        self.assertEqual(s1.shape, (_NL, _NS))

    def test_h2d_history_length(self):
        inv = _make_hybrid2d(epochs=5)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._history), 5)

    # ── HybridInverter3D ─────────────────────────────

    def test_h3d_fit_sets_is_fitted(self):
        inv = _make_hybrid3d(epochs=5)
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)

    def test_h3d_volume_shape(self):
        inv = _make_hybrid3d(epochs=5)
        inv.fit(verbose=False)
        vol = inv.resistivity_volume()
        self.assertEqual(vol.shape, (_NL, _NS))

    def test_h3d_stage1_volume_shape(self):
        inv = _make_hybrid3d(epochs=5)
        inv.fit(verbose=False)
        s1 = inv.stage1_volume()
        self.assertEqual(s1.shape, (_NL, _NS))

    def test_h3d_history_length(self):
        inv = _make_hybrid3d(epochs=4)
        inv.fit(verbose=False)
        self.assertEqual(len(inv._history), 4)


if __name__ == "__main__":
    unittest.main()
