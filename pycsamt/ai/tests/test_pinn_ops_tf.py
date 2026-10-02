# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.ai.inversion._pinn_ops_tf.

Requires TensorFlow to be installed and its native runtime to actually
load — run this file in an environment that has it (e.g. the
``base-attentiv`` conda env); it is skipped everywhere else.

TensorFlow can be "installed" (importable per ``find_spec``) while its
native DLL is broken; importing it in-process to probe availability can
trigger a native access violation that crashes the whole pytest process
outright rather than raising a catchable Python exception. The probe
below runs in a subprocess to isolate any such crash — see
``pycsamt/backends/tests/test_tensorflow_backend.py`` for the same
pattern. Do not replace this with ``pytest.importorskip("tensorflow")``.
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

_NF = 8
_NS = 3
_NL = 3
_FREQS = np.logspace(3, 0, _NF)


def _obs1d(name="S01", rho_val=100.0, phase_val=45.0):
    from pycsamt.ai.inversion._sites_bridge import SiteObs1D

    rho = np.full(_NF, rho_val)
    ph = np.full(_NF, phase_val)
    return SiteObs1D(
        name=name,
        freq=_FREQS.copy(),
        rho_obs=rho,
        phase_obs=ph,
    )


def _panel(n_stations=_NS, n_freqs=_NF, with_nan=False):
    log_rho_obs = np.full((n_stations, n_freqs), 2.0)
    ph_obs = np.full((n_stations, n_freqs), 45.0)
    if with_nan:
        log_rho_obs[0, 0] = np.nan
        ph_obs[0, 1] = np.nan
    return log_rho_obs, ph_obs


def _adjacency(n=_NS):
    A = np.ones((n, n)) - np.eye(n)
    return A.astype(np.float64)


# ── low-level physics kernels ────────────────────────────────────────────


class TestMt1dTf:
    def test_homogeneous_halfspace_rho_a_equals_true_rho(self):
        import tensorflow as tf

        from pycsamt.ai.inversion._pinn_ops_tf import _mt1d_tf

        true_rho = 100.0
        log_rho = tf.constant(
            [np.log10(true_rho), np.log10(true_rho)], dtype=tf.float64
        )
        log_thick = tf.constant([2.0], dtype=tf.float64)
        freqs_t = tf.constant(_FREQS, dtype=tf.float64)

        log10_rho_a, phase = _mt1d_tf(freqs_t, log_rho, log_thick)

        rho_a = 10.0 ** log10_rho_a.numpy()
        np.testing.assert_allclose(rho_a, true_rho, rtol=1e-6)
        np.testing.assert_allclose(phase.numpy(), 45.0, atol=1e-6)

    def test_two_layer_shapes(self):
        import tensorflow as tf

        from pycsamt.ai.inversion._pinn_ops_tf import _mt1d_tf

        log_rho = tf.constant([1.0, 2.0, 3.0], dtype=tf.float64)
        log_thick = tf.constant([1.5, 2.0], dtype=tf.float64)
        freqs_t = tf.constant(_FREQS, dtype=tf.float64)

        log10_rho_a, phase = _mt1d_tf(freqs_t, log_rho, log_thick)
        assert log10_rho_a.shape == (_NF,)
        assert phase.shape == (_NF,)
        assert np.all(np.isfinite(log10_rho_a.numpy()))
        assert np.all(np.isfinite(phase.numpy()))


class TestMt1dTfBatch:
    def test_homogeneous_halfspace_batch(self):
        import tensorflow as tf

        from pycsamt.ai.inversion._pinn_ops_tf import _mt1d_tf_batch

        true_rho = 50.0
        log_rho = tf.constant(
            np.full((_NS, _NL), np.log10(true_rho)), dtype=tf.float64
        )
        log_thick = tf.constant(
            np.full((_NS, _NL - 1), 2.0), dtype=tf.float64
        )
        freqs_t = tf.constant(_FREQS, dtype=tf.float64)

        log10_rho_a, phase = _mt1d_tf_batch(freqs_t, log_rho, log_thick)
        assert log10_rho_a.shape == (_NS, _NF)
        assert phase.shape == (_NS, _NF)

        rho_a = 10.0 ** log10_rho_a.numpy()
        np.testing.assert_allclose(rho_a, true_rho, rtol=1e-6)
        np.testing.assert_allclose(phase.numpy(), 45.0, atol=1e-6)

    def test_heterogeneous_stations_differ(self):
        import tensorflow as tf

        from pycsamt.ai.inversion._pinn_ops_tf import _mt1d_tf_batch

        log_rho = tf.constant(
            [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0], [2.0, 2.0, 2.0]],
            dtype=tf.float64,
        )
        log_thick = tf.constant(
            np.full((3, 2), 2.0), dtype=tf.float64
        )
        freqs_t = tf.constant(_FREQS, dtype=tf.float64)

        log10_rho_a, phase = _mt1d_tf_batch(freqs_t, log_rho, log_thick)
        assert not np.allclose(
            log10_rho_a.numpy()[0], log10_rho_a.numpy()[1]
        )


class TestClipGradsTf:
    def test_clips_large_gradients(self):
        import tensorflow as tf

        from pycsamt.ai.inversion._pinn_ops_tf import _clip_grads_tf

        grads = [
            tf.constant([100.0, 100.0], dtype=tf.float64),
            tf.constant([0.1], dtype=tf.float64),
        ]
        clipped = _clip_grads_tf(grads, max_norm=1.0)
        assert float(np.linalg.norm(clipped[0].numpy())) <= 1.0 + 1e-9
        np.testing.assert_allclose(clipped[1].numpy(), [0.1])

    def test_passes_through_none_gradients(self):
        import tensorflow as tf

        from pycsamt.ai.inversion._pinn_ops_tf import _clip_grads_tf

        grads = [None, tf.constant([1.0], dtype=tf.float64)]
        clipped = _clip_grads_tf(grads)
        assert clipped[0] is None
        assert clipped[1] is not None


# ── high-level dispatch through pycsamt.ai.inversion._pinn_ops ──────────


class TestFitStationTfDispatch:
    def test_returns_expected_keys_and_shapes(self):
        from pycsamt.ai.inversion._pinn_ops import fit_station

        obs = _obs1d()
        res = fit_station(
            obs,
            n_layers=_NL,
            depth_max=500.0,
            lam=0.01,
            lr=0.05,
            epochs=3,
            device="/CPU:0",
            log_every=0,
            backend="tensorflow",
        )
        assert set(res) == {"log_rho", "log_thick", "history"}
        assert res["log_rho"].shape == (_NL,)
        assert res["log_thick"].shape == (_NL - 1,)
        assert len(res["history"]) == 3
        assert all(np.isfinite(v) for v in res["history"])

    def test_verbose_logging_branch(self, capsys):
        from pycsamt.ai.inversion._pinn_ops import fit_station

        obs = _obs1d()
        fit_station(
            obs,
            n_layers=_NL,
            depth_max=500.0,
            lam=0.01,
            lr=0.05,
            epochs=2,
            device="/CPU:0",
            log_every=1,
            backend="tensorflow",
        )
        captured = capsys.readouterr()
        assert obs.name in captured.out
        assert "loss=" in captured.out

    def test_explicit_init_values_used(self):
        from pycsamt.ai.inversion._pinn_ops import fit_station

        obs = _obs1d()
        init_lr = np.array([1.0, 2.0, 3.0])
        init_lt = np.array([1.0, 1.0])
        res = fit_station(
            obs,
            n_layers=_NL,
            depth_max=500.0,
            lam=0.01,
            lr=0.0,
            epochs=1,
            device="/CPU:0",
            log_every=0,
            init_log_rho=init_lr,
            init_log_thick=init_lt,
            backend="tensorflow",
        )
        assert res["log_rho"].shape == (_NL,)

    def test_handles_missing_observations(self):
        from pycsamt.ai.inversion._sites_bridge import SiteObs1D
        from pycsamt.ai.inversion._pinn_ops import fit_station

        rho = np.full(_NF, 100.0)
        ph = np.full(_NF, 45.0)
        rho[0] = np.nan
        ph[1] = np.nan
        obs = SiteObs1D(
            name="S_nan", freq=_FREQS.copy(), rho_obs=rho, phase_obs=ph
        )
        res = fit_station(
            obs,
            n_layers=_NL,
            depth_max=500.0,
            lam=0.01,
            lr=0.05,
            epochs=2,
            device="/CPU:0",
            log_every=0,
            backend="tensorflow",
        )
        assert np.all(np.isfinite(res["log_rho"]))


    def test_all_missing_observations_falls_back_to_default_init(self):
        from pycsamt.ai.inversion._sites_bridge import SiteObs1D
        from pycsamt.ai.inversion._pinn_ops import fit_station

        rho = np.full(_NF, np.nan)
        ph = np.full(_NF, np.nan)
        obs = SiteObs1D(
            name="S_all_nan",
            freq=_FREQS.copy(),
            rho_obs=rho,
            phase_obs=ph,
        )
        res = fit_station(
            obs,
            n_layers=_NL,
            depth_max=500.0,
            lam=0.01,
            lr=0.05,
            epochs=2,
            device="/CPU:0",
            log_every=0,
            backend="tensorflow",
        )
        assert np.all(np.isfinite(res["log_rho"]))
        assert np.all(np.isfinite(res["log_thick"]))


class TestFit2dJointTfDispatch:
    def test_multi_station_shapes_and_history(self):
        from pycsamt.ai.inversion._pinn_ops import fit_2d_joint

        log_rho_obs, ph_obs = _panel()
        res = fit_2d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_x=0.01,
            lr=0.05,
            epochs=3,
            device="/CPU:0",
            log_every=0,
            backend="tensorflow",
        )
        assert res["log_rho"].shape == (_NS, _NL)
        assert res["log_thick"].shape == (_NS, _NL - 1)
        assert len(res["history"]) == 3

    def test_single_station_zero_lateral_branch(self):
        from pycsamt.ai.inversion._pinn_ops import fit_2d_joint

        log_rho_obs, ph_obs = _panel(n_stations=1)
        res = fit_2d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_x=0.5,
            lr=0.05,
            epochs=2,
            device="/CPU:0",
            log_every=0,
            backend="tensorflow",
        )
        assert res["log_rho"].shape == (1, _NL)

    def test_missing_observations_masked(self):
        from pycsamt.ai.inversion._pinn_ops import fit_2d_joint

        log_rho_obs, ph_obs = _panel(with_nan=True)
        res = fit_2d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_x=0.01,
            lr=0.05,
            epochs=2,
            device="/CPU:0",
            log_every=0,
            backend="tensorflow",
        )
        assert np.all(np.isfinite(res["log_rho"]))

    def test_verbose_progress_bar_branch(self):
        from pycsamt.ai.inversion._pinn_ops import fit_2d_joint

        log_rho_obs, ph_obs = _panel()
        res = fit_2d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_x=0.01,
            lr=0.05,
            epochs=2,
            device="/CPU:0",
            log_every=1,
            backend="tensorflow",
            verbose=0,
        )
        assert len(res["history"]) == 2

    def test_explicit_init_values(self):
        from pycsamt.ai.inversion._pinn_ops import fit_2d_joint

        log_rho_obs, ph_obs = _panel()
        init_lr = np.full((_NS, _NL), 1.5)
        init_lt = np.full((_NS, _NL - 1), 1.0)
        res = fit_2d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_x=0.01,
            lr=0.0,
            epochs=1,
            device="/CPU:0",
            log_every=0,
            init_log_rho=init_lr,
            init_log_thick=init_lt,
            backend="tensorflow",
        )
        assert res["log_rho"].shape == (_NS, _NL)


    def test_all_missing_observations_falls_back_to_default_init(self):
        from pycsamt.ai.inversion._pinn_ops import fit_2d_joint

        log_rho_obs = np.full((_NS, _NF), np.nan)
        ph_obs = np.full((_NS, _NF), np.nan)
        res = fit_2d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_x=0.01,
            lr=0.05,
            epochs=2,
            device="/CPU:0",
            log_every=0,
            backend="tensorflow",
        )
        assert np.all(np.isfinite(res["log_rho"]))


class TestFit3dJointTfDispatch:
    def test_multi_station_shapes_and_history(self):
        from pycsamt.ai.inversion._pinn_ops import fit_3d_joint

        log_rho_obs, ph_obs = _panel()
        adjacency = _adjacency()
        res = fit_3d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            adjacency,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_g=0.01,
            lr=0.05,
            epochs=3,
            device="/CPU:0",
            log_every=0,
            backend="tensorflow",
        )
        assert res["log_rho"].shape == (_NS, _NL)
        assert res["log_thick"].shape == (_NS, _NL - 1)
        assert len(res["history"]) == 3

    def test_missing_observations_masked(self):
        from pycsamt.ai.inversion._pinn_ops import fit_3d_joint

        log_rho_obs, ph_obs = _panel(with_nan=True)
        adjacency = _adjacency()
        res = fit_3d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            adjacency,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_g=0.01,
            lr=0.05,
            epochs=2,
            device="/CPU:0",
            log_every=0,
            backend="tensorflow",
        )
        assert np.all(np.isfinite(res["log_rho"]))

    def test_explicit_init_values(self):
        from pycsamt.ai.inversion._pinn_ops import fit_3d_joint

        log_rho_obs, ph_obs = _panel()
        adjacency = _adjacency()
        init_lr = np.full((_NS, _NL), 1.5)
        init_lt = np.full((_NS, _NL - 1), 1.0)
        res = fit_3d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            adjacency,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_g=0.01,
            lr=0.0,
            epochs=1,
            device="/CPU:0",
            log_every=0,
            init_log_rho=init_lr,
            init_log_thick=init_lt,
            backend="tensorflow",
        )
        assert res["log_rho"].shape == (_NS, _NL)

    def test_verbose_progress_bar_branch(self):
        from pycsamt.ai.inversion._pinn_ops import fit_3d_joint

        log_rho_obs, ph_obs = _panel()
        adjacency = _adjacency()
        res = fit_3d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            adjacency,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_g=0.01,
            lr=0.05,
            epochs=2,
            device="/CPU:0",
            log_every=1,
            backend="tensorflow",
            verbose=0,
        )
        assert len(res["history"]) == 2

    def test_all_missing_observations_falls_back_to_default_init(self):
        from pycsamt.ai.inversion._pinn_ops import fit_3d_joint

        log_rho_obs = np.full((_NS, _NF), np.nan)
        ph_obs = np.full((_NS, _NF), np.nan)
        adjacency = _adjacency()
        res = fit_3d_joint(
            log_rho_obs,
            ph_obs,
            _FREQS,
            adjacency,
            n_layers=_NL,
            depth_max=500.0,
            lam_z=0.01,
            lam_g=0.01,
            lr=0.05,
            epochs=2,
            device="/CPU:0",
            log_every=0,
            backend="tensorflow",
        )
        assert np.all(np.isfinite(res["log_rho"]))
