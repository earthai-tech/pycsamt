# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
r"""
Coverage tests for the AI-inversion Hybrid/PINN 2-D/3-D
stack and its ``_sites_bridge`` adapter layer.

Complements ``test_pinn_hybrid_inverters.py`` (which bypasses
``__init__`` via ``object.__new__`` factories) by exercising the
*real* constructors, ``_sites_bridge.py`` conversion functions
against real EDI data, and the resize / branch-selection paths
that the factory-based tests never reach.
"""

from __future__ import annotations

import unittest
import warnings
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

_HERE = Path(__file__).resolve()
_REPO_ROOT = _HERE.parents[3]
_EDI_DIR = _REPO_ROOT / "data" / "3edis"

_NF = 8
_NL = 3


def _has_torch():
    try:
        import torch  # noqa: F401

        return True
    except ImportError:
        return False


def _set_backend(name):
    from pycsamt.backends import set_backend

    set_backend(name)


_SKIP_TORCH = unittest.skipUnless(_has_torch(), "torch not available")
_SKIP_EDI = unittest.skipUnless(
    _EDI_DIR.exists(), f"sample EDI dir missing: {_EDI_DIR}"
)


# ════════════════════════════════════════════════════
# 1.  _sites_bridge.py — real-EDI conversions
# ════════════════════════════════════════════════════


@_SKIP_EDI
class TestSitesBridgeReal(unittest.TestCase):
    r"""Exercise sites_to_* against the real 3-EDI sample dir."""

    def test_sites_to_obs_1d_basic(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_1d

        obs = sites_to_obs_1d(str(_EDI_DIR), comp="xy")
        self.assertEqual(len(obs), 3)
        for o in obs:
            self.assertTrue(np.all(np.diff(o.freq) < 0))
            self.assertTrue(np.all(o.rho_obs > 0))

    def test_sites_to_obs_1d_bad_comp(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_1d

        with self.assertRaises(ValueError):
            sites_to_obs_1d(str(_EDI_DIR), comp="zz")

    def test_sites_to_obs_1d_passthrough(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_1d

        obs = sites_to_obs_1d(str(_EDI_DIR))
        obs2 = sites_to_obs_1d(obs)
        self.assertIs(obs2, obs)

    def test_sites_to_features_1d_shape(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_features_1d

        X, freqs, names = sites_to_features_1d(str(_EDI_DIR), n_freqs=16)
        self.assertEqual(X.shape, (3, 32))
        self.assertEqual(freqs.shape, (16,))
        self.assertEqual(len(names), 3)

    def test_obs_to_features_1d_fills_defaults(self):
        from pycsamt.ai.inversion._sites_bridge import (
            obs_to_features_1d,
            sites_to_obs_1d,
        )

        obs = sites_to_obs_1d(str(_EDI_DIR))
        X, freqs, names = obs_to_features_1d(obs, n_freqs=16)
        self.assertEqual(X.shape, (3, 32))
        self.assertTrue(np.all(np.isfinite(X)))

    def test_obs_to_features_1d_short_series_all_nan_filled(self):
        from pycsamt.ai.inversion._sites_bridge import (
            SiteObs1D,
            obs_to_features_1d,
        )

        obs = [
            SiteObs1D(
                name="S1",
                freq=np.array([100.0]),
                rho_obs=np.array([50.0]),
                phase_obs=np.array([30.0]),
            )
        ]
        X, freqs, names = obs_to_features_1d(obs, n_freqs=5)
        self.assertEqual(X.shape, (1, 10))
        # log-rho block filled with the 2.0 fallback constant
        self.assertTrue(np.all(X[:, :5] == 2.0))
        self.assertTrue(np.all(X[:, 5:] == 45.0))

    def test_sites_to_features_1d_single_point_series_all_nan(self):
        from pycsamt.ai.inversion._sites_bridge import (
            SiteObs1D,
            sites_to_features_1d,
        )

        obs = [
            SiteObs1D(
                name="S1",
                freq=np.array([100.0]),
                rho_obs=np.array([50.0]),
                phase_obs=np.array([30.0]),
            )
        ]
        X, freqs, names = sites_to_features_1d(obs, n_freqs=5)
        self.assertTrue(np.all(np.isnan(X)))

    def test_sites_to_obs_2d_basic(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        obs = sites_to_obs_2d(str(_EDI_DIR))
        self.assertEqual(len(obs), 3)
        for o in obs:
            self.assertEqual(o.rho_te.shape, o.rho_tm.shape)

    def test_sites_to_obs_2d_bad_comp_te(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        with self.assertRaises(ValueError):
            sites_to_obs_2d(str(_EDI_DIR), comp_te="bad")

    def test_sites_to_obs_2d_bad_comp_tm(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        with self.assertRaises(ValueError):
            sites_to_obs_2d(str(_EDI_DIR), comp_tm="bad")

    def test_sites_to_obs_2d_passthrough(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        obs = sites_to_obs_2d(str(_EDI_DIR))
        obs2 = sites_to_obs_2d(obs)
        self.assertIs(obs2, obs)

    def test_sites_to_panel_2d_shape_4ch(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_panel_2d

        panel, freqs, names = sites_to_panel_2d(str(_EDI_DIR), n_freqs=16)
        self.assertEqual(panel.shape, (1, 4, 16, 3))
        self.assertEqual(len(names), 3)

    def test_sites_to_panel_2d_shape_2ch(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_panel_2d

        panel, freqs, names = sites_to_panel_2d(
            str(_EDI_DIR), n_freqs=16, n_components=2
        )
        self.assertEqual(panel.shape[1], 2)

    def test_sites_to_panel_2d_passthrough_obs_list(self):
        from pycsamt.ai.inversion._sites_bridge import (
            sites_to_obs_2d,
            sites_to_panel_2d,
        )

        obs = sites_to_obs_2d(str(_EDI_DIR))
        panel, freqs, names = sites_to_panel_2d(obs, n_freqs=16)
        self.assertEqual(panel.shape, (1, 4, 16, 3))

    def test_sites_to_coords_3d_real(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_coords_3d

        xy = sites_to_coords_3d(str(_EDI_DIR))
        self.assertEqual(xy.shape, (3, 2))
        self.assertTrue(np.all(np.isfinite(xy)))


# ════════════════════════════════════════════════════
# 2.  _sites_bridge.py — fake-site edge cases
# ════════════════════════════════════════════════════


class _FakeSite:
    def __init__(self, name, freq=None, rho=None, phase=None, coords=None):
        self.name = name
        self.freq = freq
        self.rho = rho
        self.phase = phase
        self.coords = coords


class TestSitesBridgeFakeEdgeCases(unittest.TestCase):
    r"""
    Drive ``_extract_rho_phase``/``_extract_rho_phase_2d``/
    ``sites_to_coords_3d`` branches that need site inputs real
    EDI data cannot easily produce (bad shapes, missing arrays,
    partial coordinates).
    """

    def setUp(self):
        self._patcher = patch(
            "pycsamt.ai.inversion._sites_bridge._normalize_sites"
        )
        self._mock_normalize = self._patcher.start()
        self.addCleanup(self._patcher.stop)

    def _use(self, sites):
        self._mock_normalize.return_value = sites

    # ── 1-D extraction ───────────────────────────────

    def test_1d_missing_arrays_skipped_then_raises(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_1d

        self._use([_FakeSite("S1")])
        with self.assertRaises(ValueError):
            sites_to_obs_1d("ignored")

    def test_1d_bad_rho_ndim_warns_and_skips(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_1d

        freq = np.array([100.0, 50.0, 10.0])
        rho = np.zeros((3, 3))
        phase = np.zeros((3, 3))
        self._use([_FakeSite("S1", freq=freq, rho=rho, phase=phase)])
        with self.assertWarns(UserWarning):
            with self.assertRaises(ValueError):
                sites_to_obs_1d("ignored")

    def test_1d_rho_1d_path_and_sorting(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_1d

        freq = np.array([10.0, 100.0, 50.0])
        rho = np.array([30.0, 10.0, 20.0])
        phase = np.array([5.0, 45.0, 25.0])
        self._use([_FakeSite("S1", freq=freq, rho=rho, phase=phase)])
        out = sites_to_obs_1d("ignored")
        self.assertEqual(len(out), 1)
        o = out[0]
        np.testing.assert_allclose(o.freq, [100.0, 50.0, 10.0])
        np.testing.assert_allclose(o.rho_obs, [10.0, 20.0, 30.0])

    def test_1d_invalid_values_filtered(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_1d

        freq = np.array([100.0, 50.0, 10.0])
        rho = np.array([10.0, -5.0, np.nan])
        phase = np.array([45.0, 20.0, 5.0])
        self._use([_FakeSite("S1", freq=freq, rho=rho, phase=phase)])
        out = sites_to_obs_1d("ignored")
        self.assertEqual(len(out[0].freq), 1)

    def test_1d_all_values_invalid_skipped_silently(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_1d

        good = _FakeSite(
            "S1",
            freq=np.array([100.0, 10.0]),
            rho=np.array([10.0, 20.0]),
            phase=np.array([30.0, 40.0]),
        )
        all_bad = _FakeSite(
            "S2",
            freq=np.array([100.0, 10.0]),
            rho=np.array([-1.0, np.nan]),
            phase=np.array([30.0, 40.0]),
        )
        self._use([good, all_bad])
        out = sites_to_obs_1d("ignored", verbose=0)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0].name, "S1")

    def test_1d_verbose_warns_on_partial_skip(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_1d

        good = _FakeSite(
            "S1",
            freq=np.array([100.0, 10.0]),
            rho=np.array([10.0, 20.0]),
            phase=np.array([30.0, 40.0]),
        )
        bad = _FakeSite("S2")
        self._use([good, bad])
        with self.assertWarns(UserWarning):
            out = sites_to_obs_1d("ignored", verbose=1)
        self.assertEqual(len(out), 1)

    # ── obs_to_features_1d ────────────────────────────

    def test_obs_to_features_1d_skips_objects_without_rho(self):
        from pycsamt.ai.inversion._sites_bridge import (
            SiteObs1D,
            obs_to_features_1d,
        )

        class _NoRho:
            def __init__(self, name, freq):
                self.name = name
                self.freq = freq

        good = SiteObs1D(
            name="S1",
            freq=np.array([100.0, 10.0]),
            rho_obs=np.array([10.0, 20.0]),
            phase_obs=np.array([30.0, 40.0]),
        )
        bad = _NoRho("S2", np.array([100.0, 10.0]))
        X, freqs, names = obs_to_features_1d([good, bad], n_freqs=4)
        self.assertEqual(names, ["S1"])
        self.assertEqual(X.shape, (1, 8))

    def test_obs_to_features_1d_all_skipped_raises(self):
        from pycsamt.ai.inversion._sites_bridge import obs_to_features_1d

        class _NoRho:
            def __init__(self, name, freq):
                self.name = name
                self.freq = freq

        with self.assertRaises(ValueError):
            obs_to_features_1d(
                [_NoRho("S1", np.array([100.0, 10.0]))], n_freqs=4
            )

    # ── 2-D extraction ───────────────────────────────

    def test_2d_missing_arrays_skipped_then_raises(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        self._use([_FakeSite("S1")])
        with self.assertRaises(ValueError):
            sites_to_obs_2d("ignored")

    def test_2d_te_and_tm_both_invalid_skipped(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        n = 2
        freq = np.array([100.0, 50.0])
        rho = np.ones((n, 2, 2)) * -1.0  # every component invalid
        phase = np.ones((n, 2, 2)) * 30.0
        self._use([_FakeSite("S1", freq=freq, rho=rho, phase=phase)])
        with self.assertRaises(ValueError):
            sites_to_obs_2d("ignored")

    def test_2d_verbose_warns_on_partial_skip(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        n = 2
        freq = np.array([100.0, 50.0])
        good = _FakeSite(
            "S1",
            freq=freq,
            rho=np.ones((n, 2, 2)) * 10.0,
            phase=np.ones((n, 2, 2)) * 30.0,
        )
        bad = _FakeSite(
            "S2",
            freq=freq,
            rho=np.ones((n, 2, 2)) * -1.0,
            phase=np.ones((n, 2, 2)) * 30.0,
        )
        self._use([good, bad])
        with self.assertWarns(UserWarning):
            out = sites_to_obs_2d("ignored", verbose=1)
        self.assertEqual(len(out), 1)

    def test_2d_rho_1d_path_tm_mirrors_te(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        freq = np.array([100.0, 50.0])
        rho = np.array([10.0, 20.0])
        phase = np.array([-30.0, -40.0])
        self._use([_FakeSite("S1", freq=freq, rho=rho, phase=phase)])
        out = sites_to_obs_2d("ignored")
        o = out[0]
        np.testing.assert_allclose(o.rho_tm, o.rho_te)
        np.testing.assert_allclose(o.phase_tm, np.abs(o.phase_te))

    def test_2d_bad_rho_ndim_skipped_no_warning(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        freq = np.array([100.0, 50.0])
        rho = np.zeros((2, 3))
        phase = np.zeros((2, 3))
        good = _FakeSite(
            "S1",
            freq=freq,
            rho=np.ones((2, 2, 2)),
            phase=np.ones((2, 2, 2)) * 45.0,
        )
        bad = _FakeSite("S2", freq=freq, rho=rho, phase=phase)
        self._use([good, bad])
        out = sites_to_obs_2d("ignored")
        self.assertEqual(len(out), 1)

    def test_2d_tm_invalid_falls_back_to_te(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        n = 4
        freq = np.array([100.0, 50.0, 10.0, 5.0])
        rho = np.ones((n, 2, 2)) * 10.0
        phase = np.ones((n, 2, 2)) * 30.0
        # yx (TM, comp index (1,0)) invalid on 2 of 4 rows
        rho[1, 1, 0] = -1.0
        rho[2, 1, 0] = np.nan
        self._use([_FakeSite("S1", freq=freq, rho=rho, phase=phase)])
        out = sites_to_obs_2d("ignored")
        o = out[0]
        # TE-valid mask keeps all 4 rows; TM filled from TE where invalid
        self.assertEqual(len(o.freq), 4)
        np.testing.assert_allclose(o.rho_tm[1], o.rho_te[1])
        np.testing.assert_allclose(o.rho_tm[2], o.rho_te[2])

    def test_2d_te_invalid_uses_tm_mask(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_obs_2d

        n = 3
        freq = np.array([100.0, 50.0, 10.0])
        rho = np.ones((n, 2, 2)) * 10.0
        phase = np.ones((n, 2, 2)) * 30.0
        # xy (TE) invalid everywhere; yx (TM) stays valid
        rho[:, 0, 1] = -1.0
        self._use([_FakeSite("S1", freq=freq, rho=rho, phase=phase)])
        out = sites_to_obs_2d("ignored")
        self.assertEqual(len(out[0].freq), 3)

    def test_sites_to_panel_2d_no_valid_raises(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_panel_2d

        self._use([_FakeSite("S1")])
        with self.assertRaises(ValueError):
            sites_to_panel_2d("ignored")

    # ── coords_3d ─────────────────────────────────────

    def test_coords_3d_fallback_uniform_grid(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_coords_3d

        self._use([_FakeSite(f"S{i}") for i in range(4)])
        xy = sites_to_coords_3d("ignored", station_spacing=100.0)
        expected = np.array(
            [[0.0, 0.0], [100.0, 0.0], [0.0, 100.0], [100.0, 100.0]]
        )
        np.testing.assert_allclose(xy, expected)

    def test_coords_3d_bad_coords_type_caught(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_coords_3d

        self._use(
            [
                _FakeSite("S1", coords=("bad", "data")),
                _FakeSite("S2", coords=(1.0, 2.0)),
            ]
        )
        xy = sites_to_coords_3d("ignored")
        self.assertEqual(xy.shape, (2, 2))
        self.assertTrue(np.all(np.isfinite(xy)))

    def test_coords_3d_partial_missing_filled_zero(self):
        from pycsamt.ai.inversion._sites_bridge import sites_to_coords_3d

        sites = [
            _FakeSite("S1", coords=(10.0, 10.0)),
            _FakeSite("S2", coords=(10.01, 10.01)),
            _FakeSite("S3", coords=(10.02, 10.0)),
            _FakeSite("S4", coords=None),
        ]
        self._use(sites)
        xy = sites_to_coords_3d("ignored")
        self.assertEqual(xy.shape, (4, 2))
        # last site had no coords -> filled with 0.0 in the flat-earth branch
        np.testing.assert_allclose(xy[3], [0.0, 0.0])
        self.assertTrue(np.all(np.isfinite(xy)))


# ════════════════════════════════════════════════════
# 3.  Real __init__ + fit() for Hybrid2D/3D, PINN2D/3D
# ════════════════════════════════════════════════════


def _fake_ai2d(n_sta, n_depth=_NL, n_ch=2, n_f=_NF):
    from pycsamt.ai.inversion.inv2d import EMInverter2D

    ai = EMInverter2D(
        n_components=n_ch, n_depth=n_depth, n_stations=n_sta, n_freqs=n_f
    )
    ai._is_fitted = True

    def _predict(panel, as_log_rho=False):
        n_sta_out = panel.shape[3]
        return np.full((1, n_depth, n_sta_out), 2.0)

    ai.predict = _predict
    return ai


def _synth_obs2d(n=3, nf=_NF):
    r"""
    Same-frequency-grid synthetic SiteObs2D list.

    Real EDI stations have differing frequency ranges, so the
    auto-detected common grid extrapolates to NaN at the edges
    for at least one station. ``fit_2d_joint``/``fit_3d_joint``
    (``_pinn_ops_torch.py``, out of scope for this batch) do not
    mask that NaN cleanly and the whole optimisation comes back
    NaN -- a real, but out-of-scope, bug. Fit-path tests use this
    NaN-free synthetic panel instead so they exercise the Hybrid/
    PINN classes' own logic without tripping that landmine.
    """
    from pycsamt.ai.inversion._sites_bridge import SiteObs2D

    freqs = np.logspace(3, 0, nf)
    out = []
    for i in range(n):
        rho = np.full(nf, 50.0 + 10.0 * i)
        ph = np.full(nf, 45.0)
        out.append(
            SiteObs2D(
                name=f"S{i + 1}",
                freq=freqs.copy(),
                rho_te=rho.copy(),
                phase_te=ph.copy(),
                rho_tm=rho.copy() * 1.1,
                phase_tm=ph.copy(),
            )
        )
    return out


def _fake_ai3d(n_sta, n_layers=_NL, n_features=None):
    from pycsamt.ai.inversion.inv3d import GCNInverter3D

    nf = n_features if n_features is not None else _NF * 2
    ai = GCNInverter3D(n_features=nf, n_layers=n_layers)
    ai._is_fitted = True

    def _predict(X, adjacency=None, as_log_rho=False):
        n = X.shape[0]
        return np.full((n, 2 * n_layers - 1), 2.0)

    ai.predict = _predict
    return ai


@_SKIP_TORCH
@_SKIP_EDI
class TestHybrid2DRealConstruction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _set_backend("torch")

    @classmethod
    def tearDownClass(cls):
        _set_backend("auto")

    def test_mode_invalid_raises(self):
        from pycsamt.ai.inversion import HybridInverter2D

        with self.assertRaises(ValueError):
            HybridInverter2D(
                str(_EDI_DIR), _fake_ai2d(3), mode="bad", n_freqs=_NF
            )

    def test_ai_inverter_bad_type_raises(self):
        from pycsamt.ai.inversion import HybridInverter2D

        with self.assertRaises(TypeError):
            HybridInverter2D(str(_EDI_DIR), object(), n_freqs=_NF)

    def test_ai_inverter_unfitted_raises(self):
        from pycsamt.ai.inversion import HybridInverter2D

        ai = _fake_ai2d(3)
        ai._is_fitted = False
        with self.assertRaises(ValueError):
            HybridInverter2D(str(_EDI_DIR), ai, n_freqs=_NF)

    def test_ai_inverter_loaded_from_path(self):
        from pycsamt.ai.inversion import HybridInverter2D
        from pycsamt.ai.inversion.inv2d import EMInverter2D

        fake = _fake_ai2d(3)
        with patch.object(EMInverter2D, "load", return_value=fake):
            inv = HybridInverter2D(
                str(_EDI_DIR), "checkpoint.npz", n_freqs=_NF
            )
        self.assertIs(inv._ai_inv, fake)

    def test_init_matching_station_count(self):
        from pycsamt.ai.inversion import HybridInverter2D

        inv = HybridInverter2D(
            str(_EDI_DIR), _fake_ai2d(3), n_freqs=_NF, epochs=3
        )
        self.assertEqual(inv.n_sites, 3)
        self.assertEqual(len(inv.stations), 3)
        self.assertIn("unfitted", repr(inv))

    def test_fit_with_resize_mismatched_ai_stations(self):
        from pycsamt.ai.inversion import HybridInverter2D

        # ai trained for 5 stations, obs has 3 -> both panel and
        # section resize helpers get exercised.
        inv = HybridInverter2D(
            _synth_obs2d(3),
            _fake_ai2d(5),
            n_freqs=_NF,
            epochs=3,
            mode="both",
        )
        inv.fit(verbose=False)
        self.assertTrue(inv._is_fitted)
        self.assertIn("fitted", repr(inv))
        sec = inv.resistivity_section(as_log10=False)
        self.assertEqual(sec.shape, (_NL, 3))
        self.assertTrue(np.all(sec > 0))
        th = inv.thickness_section()
        self.assertEqual(th.shape, (_NL - 1, 3))
        s1 = inv.stage1_section(as_log10=False)
        self.assertEqual(s1.shape, (_NL, 3))

    def test_stage1_section_before_fit_raises(self):
        from pycsamt.ai.inversion import HybridInverter2D

        inv = HybridInverter2D(
            str(_EDI_DIR), _fake_ai2d(3), n_freqs=_NF, epochs=3
        )
        with self.assertRaises(RuntimeError):
            inv.stage1_section()

    def test_residuals_stage1_and_stage2_tm_mode(self):
        from pycsamt.ai.inversion import HybridInverter2D

        inv = HybridInverter2D(
            _synth_obs2d(3),
            _fake_ai2d(3),
            n_freqs=_NF,
            epochs=3,
            mode="tm",
        )
        inv.fit(verbose=False)
        df1 = inv.residuals(stage=1)
        df2 = inv.residuals(stage=2)
        self.assertGreater(len(df1), 0)
        self.assertGreater(len(df2), 0)
        self.assertIn("rho_obs", df1.columns)

    def test_convergence_curve_df(self):
        import pandas as pd

        from pycsamt.ai.inversion import HybridInverter2D

        inv = HybridInverter2D(
            _synth_obs2d(3), _fake_ai2d(3), n_freqs=_NF, epochs=4
        )
        inv.fit(verbose=False)
        df = inv.convergence_curve()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 4)


@_SKIP_TORCH
@_SKIP_EDI
class TestHybrid3DRealConstruction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _set_backend("torch")

    @classmethod
    def tearDownClass(cls):
        _set_backend("auto")

    def test_mode_invalid_raises(self):
        from pycsamt.ai.inversion import HybridInverter3D

        with self.assertRaises(ValueError):
            HybridInverter3D(
                str(_EDI_DIR), _fake_ai3d(3), mode="bad", n_freqs=_NF
            )

    def test_ai_inverter_bad_type_raises(self):
        from pycsamt.ai.inversion import HybridInverter3D

        with self.assertRaises(TypeError):
            HybridInverter3D(str(_EDI_DIR), object(), n_freqs=_NF)

    def test_ai_inverter_unfitted_raises(self):
        from pycsamt.ai.inversion import HybridInverter3D

        ai = _fake_ai3d(3)
        ai._is_fitted = False
        with self.assertRaises(ValueError):
            HybridInverter3D(str(_EDI_DIR), ai, n_freqs=_NF)

    def test_ai_inverter_loaded_from_path(self):
        from pycsamt.ai.inversion import HybridInverter3D
        from pycsamt.ai.inversion.inv3d import GCNInverter3D

        fake = _fake_ai3d(3)
        with patch.object(GCNInverter3D, "load", return_value=fake):
            inv = HybridInverter3D(
                str(_EDI_DIR), "checkpoint.npz", n_freqs=_NF
            )
        self.assertIs(inv._ai_inv, fake)

    def test_init_default_coords_and_adjacency(self):
        from pycsamt.ai.inversion import HybridInverter3D

        inv = HybridInverter3D(
            str(_EDI_DIR), _fake_ai3d(3), n_freqs=_NF, epochs=3
        )
        self.assertEqual(inv.station_coords().shape, (3, 2))
        self.assertEqual(inv.adjacency().shape, (3, 3))
        self.assertIn("unfitted", repr(inv))

    def test_init_explicit_coords_and_adjacency(self):
        from pycsamt.ai.inversion import HybridInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        adj = np.eye(3)
        inv = HybridInverter3D(
            str(_EDI_DIR),
            _fake_ai3d(3),
            n_freqs=_NF,
            epochs=3,
            station_coords=coords,
            adjacency=adj,
        )
        np.testing.assert_allclose(inv.station_coords(), coords)
        np.testing.assert_allclose(inv.adjacency(), adj)

    def test_fit_and_outputs_both_mode(self):
        from pycsamt.ai.inversion import HybridInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        inv = HybridInverter3D(
            _synth_obs2d(3),
            _fake_ai3d(3),
            n_freqs=_NF,
            epochs=3,
            mode="both",
            station_coords=coords,
        )
        inv.fit(verbose=False)
        self.assertIn("fitted", repr(inv))
        vol = inv.resistivity_volume(as_log10=False)
        self.assertEqual(vol.shape, (_NL, 3))
        self.assertTrue(np.all(vol > 0))
        th = inv.thickness_volume()
        self.assertEqual(th.shape, (_NL - 1, 3))
        s1 = inv.stage1_volume(as_log10=False)
        self.assertEqual(s1.shape, (_NL, 3))

    def test_stage1_volume_before_fit_raises(self):
        from pycsamt.ai.inversion import HybridInverter3D

        inv = HybridInverter3D(
            str(_EDI_DIR), _fake_ai3d(3), n_freqs=_NF, epochs=3
        )
        with self.assertRaises(RuntimeError):
            inv.stage1_volume()

    def test_residuals_stage1_and_stage2_tm_mode(self):
        from pycsamt.ai.inversion import HybridInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        inv = HybridInverter3D(
            _synth_obs2d(3),
            _fake_ai3d(3),
            n_freqs=_NF,
            epochs=3,
            mode="tm",
            station_coords=coords,
        )
        inv.fit(verbose=False)
        df1 = inv.residuals(stage=1)
        df2 = inv.residuals(stage=2)
        self.assertGreater(len(df1), 0)
        self.assertGreater(len(df2), 0)


@_SKIP_TORCH
@_SKIP_EDI
class TestPINN2DRealConstruction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _set_backend("torch")

    @classmethod
    def tearDownClass(cls):
        _set_backend("auto")

    def test_mode_invalid_raises(self):
        from pycsamt.ai.inversion import PINNInverter2D

        with self.assertRaises(ValueError):
            PINNInverter2D(str(_EDI_DIR), mode="bad", n_freqs=_NF)

    def test_fit_both_mode_outputs(self):
        from pycsamt.ai.inversion import PINNInverter2D

        inv = PINNInverter2D(
            _synth_obs2d(3),
            n_layers=_NL,
            n_freqs=_NF,
            epochs=3,
            mode="both",
        )
        inv.fit(verbose=False)
        self.assertIn("fitted", repr(inv))
        sec = inv.resistivity_section(as_log10=False)
        self.assertEqual(sec.shape, (_NL, 3))
        self.assertTrue(np.all(sec > 0))
        th = inv.thickness_section()
        self.assertEqual(th.shape, (_NL - 1, 3))

    def test_fit_tm_mode_residuals(self):
        from pycsamt.ai.inversion import PINNInverter2D

        inv = PINNInverter2D(
            _synth_obs2d(3),
            n_layers=_NL,
            n_freqs=_NF,
            epochs=3,
            mode="tm",
        )
        inv.fit(verbose=False)
        df = inv.residuals()
        self.assertGreater(len(df), 0)

    def test_init_unfitted_repr(self):
        from pycsamt.ai.inversion import PINNInverter2D

        inv = PINNInverter2D(str(_EDI_DIR), n_freqs=_NF)
        self.assertIn("unfitted", repr(inv))
        self.assertEqual(inv.n_sites, 3)


@_SKIP_TORCH
@_SKIP_EDI
class TestPINN3DRealConstruction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _set_backend("torch")

    @classmethod
    def tearDownClass(cls):
        _set_backend("auto")

    def test_mode_invalid_raises(self):
        from pycsamt.ai.inversion import PINNInverter3D

        with self.assertRaises(ValueError):
            PINNInverter3D(str(_EDI_DIR), mode="bad", n_freqs=_NF)

    def test_init_default_coords_and_adjacency(self):
        from pycsamt.ai.inversion import PINNInverter3D

        inv = PINNInverter3D(str(_EDI_DIR), n_freqs=_NF)
        self.assertEqual(inv.station_coords().shape, (3, 2))
        self.assertEqual(inv.adjacency().shape, (3, 3))

    def test_init_explicit_coords_and_adjacency(self):
        from pycsamt.ai.inversion import PINNInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        adj = np.eye(3)
        inv = PINNInverter3D(
            str(_EDI_DIR),
            n_freqs=_NF,
            station_coords=coords,
            adjacency=adj,
        )
        np.testing.assert_allclose(inv.station_coords(), coords)
        np.testing.assert_allclose(inv.adjacency(), adj)

    def test_fit_both_mode_outputs(self):
        from pycsamt.ai.inversion import PINNInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        inv = PINNInverter3D(
            _synth_obs2d(3),
            n_layers=_NL,
            n_freqs=_NF,
            epochs=3,
            mode="both",
            station_coords=coords,
        )
        inv.fit(verbose=False)
        self.assertIn("fitted", repr(inv))
        vol = inv.resistivity_volume(as_log10=False)
        self.assertEqual(vol.shape, (_NL, 3))
        th = inv.thickness_volume()
        self.assertEqual(th.shape, (_NL - 1, 3))

    def test_fit_tm_mode_residuals(self):
        from pycsamt.ai.inversion import PINNInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        inv = PINNInverter3D(
            _synth_obs2d(3),
            n_layers=_NL,
            n_freqs=_NF,
            epochs=3,
            mode="tm",
            station_coords=coords,
        )
        inv.fit(verbose=False)
        df = inv.residuals()
        self.assertGreater(len(df), 0)


# ════════════════════════════════════════════════════
# 4.  Remaining branch coverage: verbose prints,
#     residuals() exception fallback, stations
#     property, resize/padding/fill-NaN internals.
# ════════════════════════════════════════════════════


@_SKIP_TORCH
class TestRemainingBranches(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _set_backend("torch")

    @classmethod
    def tearDownClass(cls):
        _set_backend("auto")

    # ── verbose=True print branches ──────────────────

    def test_hybrid2d_fit_verbose_prints(self):
        from pycsamt.ai.inversion import HybridInverter2D

        inv = HybridInverter2D(
            _synth_obs2d(3), _fake_ai2d(3), n_freqs=_NF, epochs=2
        )
        inv.fit(verbose=True, log_every=1)
        self.assertTrue(inv._is_fitted)

    def test_hybrid3d_fit_verbose_prints(self):
        from pycsamt.ai.inversion import HybridInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        inv = HybridInverter3D(
            _synth_obs2d(3),
            _fake_ai3d(3),
            n_freqs=_NF,
            epochs=2,
            station_coords=coords,
        )
        inv.fit(verbose=True, log_every=1)
        self.assertTrue(inv._is_fitted)

    def test_pinn2d_fit_verbose_prints(self):
        from pycsamt.ai.inversion import PINNInverter2D

        inv = PINNInverter2D(
            _synth_obs2d(3), n_layers=_NL, n_freqs=_NF, epochs=2
        )
        inv.fit(verbose=True, log_every=1)
        self.assertTrue(inv._is_fitted)

    def test_pinn3d_fit_verbose_prints(self):
        from pycsamt.ai.inversion import PINNInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        inv = PINNInverter3D(
            _synth_obs2d(3),
            n_layers=_NL,
            n_freqs=_NF,
            epochs=2,
            station_coords=coords,
        )
        inv.fit(verbose=True, log_every=1)
        self.assertTrue(inv._is_fitted)

    # ── residuals() exception fallback (NaN row) ─────

    def test_hybrid2d_residuals_exception_fallback(self):
        from pycsamt.ai.inversion import HybridInverter2D

        inv = HybridInverter2D(
            _synth_obs2d(3), _fake_ai2d(3), n_freqs=_NF, epochs=2
        )
        inv.fit(verbose=False)
        with patch(
            "pycsamt.forward.em1d.MT1DForward.run",
            side_effect=RuntimeError("boom"),
        ):
            df = inv.residuals(stage=2)
        self.assertTrue(df["rho_pred"].isna().all())

    def test_hybrid3d_residuals_exception_fallback(self):
        from pycsamt.ai.inversion import HybridInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        inv = HybridInverter3D(
            _synth_obs2d(3),
            _fake_ai3d(3),
            n_freqs=_NF,
            epochs=2,
            station_coords=coords,
        )
        inv.fit(verbose=False)
        with patch(
            "pycsamt.forward.em1d.MT1DForward.run",
            side_effect=RuntimeError("boom"),
        ):
            df = inv.residuals(stage=2)
        self.assertTrue(df["rho_pred"].isna().all())

    def test_pinn2d_residuals_exception_fallback(self):
        from pycsamt.ai.inversion import PINNInverter2D

        inv = PINNInverter2D(
            _synth_obs2d(3), n_layers=_NL, n_freqs=_NF, epochs=2
        )
        inv.fit(verbose=False)
        with patch(
            "pycsamt.forward.em1d.MT1DForward.run",
            side_effect=RuntimeError("boom"),
        ):
            df = inv.residuals()
        self.assertTrue(df["rho_pred"].isna().all())

    def test_pinn3d_residuals_exception_fallback(self):
        from pycsamt.ai.inversion import PINNInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        inv = PINNInverter3D(
            _synth_obs2d(3),
            n_layers=_NL,
            n_freqs=_NF,
            epochs=2,
            station_coords=coords,
        )
        inv.fit(verbose=False)
        with patch(
            "pycsamt.forward.em1d.MT1DForward.run",
            side_effect=RuntimeError("boom"),
        ):
            df = inv.residuals()
        self.assertTrue(df["rho_pred"].isna().all())

    # ── stations property (3-D variants) ─────────────

    def test_hybrid3d_stations_property(self):
        from pycsamt.ai.inversion import HybridInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        inv = HybridInverter3D(
            _synth_obs2d(3),
            _fake_ai3d(3),
            n_freqs=_NF,
            station_coords=coords,
        )
        self.assertEqual(inv.stations, ["S1", "S2", "S3"])

    def test_pinn3d_stations_property(self):
        from pycsamt.ai.inversion import PINNInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        inv = PINNInverter3D(
            _synth_obs2d(3), n_freqs=_NF, station_coords=coords
        )
        self.assertEqual(inv.stations, ["S1", "S2", "S3"])

    # ── hybrid2d: stage-1 depth padding branch ───────

    def test_hybrid2d_stage1_padding_when_ai_depth_lt_n_layers(self):
        from pycsamt.ai.inversion import HybridInverter2D

        # ai_inverter's n_depth (2) < requested n_layers (5) ->
        # deeper layers padded with the deepest predicted value.
        inv = HybridInverter2D(
            _synth_obs2d(3),
            _fake_ai2d(3, n_depth=2),
            n_layers=5,
            n_freqs=_NF,
            epochs=2,
        )
        inv.fit(verbose=False)
        s1 = inv.stage1_section(as_log10=True)
        self.assertEqual(s1.shape, (5, 3))
        np.testing.assert_allclose(s1[3], s1[1])
        np.testing.assert_allclose(s1[4], s1[1])

    # ── hybrid2d: resize helpers direct shortcut ─────

    def test_resize_panel_stations_noop_when_equal(self):
        from pycsamt.ai.inversion.hybrid2d import _resize_panel_stations

        panel = np.random.rand(1, 2, 4, 3)
        out = _resize_panel_stations(panel, 3)
        self.assertIs(out, panel)

    def test_resize_section_stations_noop_when_equal(self):
        from pycsamt.ai.inversion.hybrid2d import _resize_section_stations

        section = np.random.rand(4, 3)
        out = _resize_section_stations(section, 3)
        self.assertIs(out, section)

    def test_fill_nan_panel_all_nan_column_fallback_zero(self):
        from pycsamt.ai.inversion.hybrid2d import _fill_nan_panel

        panel = np.full((1, 1, 4, 1), np.nan)
        out = _fill_nan_panel(panel)
        self.assertTrue(np.all(out == 0.0))

    # ── hybrid3d: uniform-thickness fallback branch ──

    def test_hybrid3d_uniform_thickness_fallback_on_shape_mismatch(self):
        from pycsamt.ai.inversion import HybridInverter3D

        coords = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
        # ai trained with n_layers=3 (outputs width 2*3-1=5) but the
        # Hybrid requests n_layers=4 -> init_log_thick width (1) !=
        # n_l-1 (3), forcing the uniform-thickness fallback.
        inv = HybridInverter3D(
            _synth_obs2d(3),
            _fake_ai3d(3, n_layers=3),
            n_layers=4,
            n_freqs=_NF,
            epochs=2,
            station_coords=coords,
        )
        inv.fit(verbose=False)
        th = inv.thickness_volume()
        self.assertEqual(th.shape, (3, 3))
        self.assertTrue(np.all(np.isfinite(th)))


if __name__ == "__main__":
    unittest.main()
