# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Coverage tests for pycsamt.ai._base.

BaseEMNet/BaseEMProcessor are abstract; concrete estimators elsewhere
(EMInverter1D, EMDenoiser, ...) exercise their ``fit``/``predict`` paths
but never touch the shared helpers exhaustively: ``_coerce_input`` (the
ndarray / Z / list-of-Z / ForwardResponse dispatch), ``_compute_metric``,
``fit_predict``, ``from_pretrained``/``load_pretrained``, and the default
save/load plumbing.  This file drives those directly through minimal
concrete subclasses.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai._base import (
    BaseEMNet,
    BaseEMProcessor,
    EMCheckpoint,
    _compute_metric,
    _z_list_to_array,
)


# ─────────────────────────────────────────────────────────────────────────────
# _compute_metric
# ─────────────────────────────────────────────────────────────────────────────


def test_compute_metric_rmse_mae_r2_relative_rmse():
    y_true = np.array([1.0, 2.0, 3.0, 4.0])
    y_pred = np.array([1.1, 1.9, 3.2, 3.8])

    rmse = _compute_metric(y_true, y_pred, "rmse")
    mae = _compute_metric(y_true, y_pred, "mae")
    r2 = _compute_metric(y_true, y_pred, "r2")
    rel = _compute_metric(y_true, y_pred, "relative_rmse")

    assert rmse == pytest.approx(
        float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    )
    assert mae == pytest.approx(float(np.mean(np.abs(y_true - y_pred))))
    assert r2 <= 1.0
    assert rel > 0.0

    # metric name is case-insensitive
    assert _compute_metric(y_true, y_pred, "RMSE") == rmse


def test_compute_metric_unknown_raises():
    with pytest.raises(ValueError, match="Unknown metric"):
        _compute_metric(np.zeros(3), np.zeros(3), "bogus")


def test_compute_metric_all_nan_returns_nan():
    y_true = np.array([np.nan, np.nan])
    y_pred = np.array([1.0, 2.0])
    assert np.isnan(_compute_metric(y_true, y_pred, "rmse"))


def test_compute_metric_ignores_nan_entries():
    y_true = np.array([1.0, np.nan, 3.0])
    y_pred = np.array([1.0, 99.0, 3.0])
    assert _compute_metric(y_true, y_pred, "rmse") == pytest.approx(0.0)


# ─────────────────────────────────────────────────────────────────────────────
# _z_list_to_array / _coerce_input — real Z objects
# ─────────────────────────────────────────────────────────────────────────────


def _make_real_z(n=6, seed=0):
    from pycsamt.z.z import Z

    rng = np.random.default_rng(seed)
    freq = np.logspace(0, 2, n)
    z_arr = np.zeros((n, 2, 2), dtype=complex)
    z_arr[:, 0, 1] = 1.0 + 0.1 * rng.standard_normal(n) + 1j * (
        1.0 + 0.1 * rng.standard_normal(n)
    )
    z_arr[:, 1, 0] = -z_arr[:, 0, 1]
    return Z(z_array=z_arr, freq=freq)


def test_z_list_to_array_real_z_single_and_list():
    z = _make_real_z()
    out_single = _z_list_to_array([z])
    assert out_single.shape[0] == 1
    assert out_single.dtype == np.float32

    z2 = _make_real_z(seed=1)
    out_list = _z_list_to_array([z, z2])
    assert out_list.shape[0] == 2
    # include_phase=True doubles the per-site width vs include_phase=False
    out_no_phase = _z_list_to_array([z, z2], include_phase=False)
    assert out_no_phase.shape[1] * 2 == out_list.shape[1]


def test_z_list_to_array_log_rho_false_skips_log10():
    z = _make_real_z()
    out_log = _z_list_to_array([z], log_rho=True)
    out_lin = _z_list_to_array([z], log_rho=False)
    n = out_lin.shape[1] // 2
    np.testing.assert_allclose(
        10.0 ** out_log[:, :n], out_lin[:, :n], rtol=1e-4
    )


def test_z_list_to_array_uncomputed_z_raises_value_error():
    from pycsamt.z.z import Z

    z = Z()  # no data attached -> resistivity/phase never computed
    with pytest.raises(ValueError, match="compute_resistivity_phase"):
        _z_list_to_array([z])


def test_z_list_to_array_duck_typed_stand_in():
    class _FakeZ:
        def __init__(self, rho_xy, phs_xy):
            self.resistivity_xy = rho_xy
            self.phase_xy = phs_xy

    fake = _FakeZ(np.array([10.0, 20.0, 30.0]), np.array([45.0, 40.0, 35.0]))
    out = _z_list_to_array([fake])
    assert out.shape == (1, 6)


# ─────────────────────────────────────────────────────────────────────────────
# BaseEMNet — minimal concrete subclass
# ─────────────────────────────────────────────────────────────────────────────


class _ToyNet(BaseEMNet):
    """A trivial linear-regression-as-a-network stand-in."""

    def _build_network(self):
        return {"w": np.zeros(self.n_layers, dtype=np.float32)}

    def fit(self, X, y=None, **kwargs):
        X = self._coerce_input(X)
        y = np.asarray(y, dtype=np.float32)
        self._network = self._build_network()
        # closed-form per-output mean (enough to make predict deterministic)
        self._network["w"] = y.mean(axis=0)
        self._is_fitted = True
        self._history = {"train_loss": [1.0, 0.5]}
        return self

    def predict(self, X):
        X = self._coerce_input(X)
        n = X.shape[0]
        return np.tile(self._network["w"], (n, 1))

    def _get_weights(self):
        return {"w": self._network["w"]}

    def _load_weights(self, weights):
        self._network = {"w": weights["w"]}


def test_base_emnet_fit_predict_score_and_fit_predict():
    net = _ToyNet(n_layers=3)
    X = np.zeros((5, 4), dtype=np.float32)
    y = np.ones((5, 3), dtype=np.float32) * 2.0

    net.fit(X, y)
    assert net._is_fitted

    score = net.score(X, y, metric="mae")
    assert score == pytest.approx(0.0, abs=1e-6)

    net2 = _ToyNet(n_layers=3)
    y_pred = net2.fit_predict(X, y)
    np.testing.assert_allclose(y_pred, y)


def test_base_emnet_save_load_round_trip(tmp_path):
    net = _ToyNet(n_layers=2, arch="toy", solver="mt1d", device="cpu")
    X = np.zeros((4, 2), dtype=np.float32)
    y = np.array([[1.0, 2.0]] * 4, dtype=np.float32)
    net.fit(X, y)

    path = tmp_path / "toy.npz"
    net.save(path)
    loaded = _ToyNet.load(path)

    assert loaded.arch == "toy"
    assert loaded.n_layers == 2
    assert loaded._is_fitted
    np.testing.assert_allclose(loaded.predict(X), net.predict(X))


def test_base_emnet_from_pretrained_not_implemented():
    with pytest.raises(NotImplementedError):
        _ToyNet.from_pretrained("mt1d-toy-v1")


def test_base_emnet_repr_fitted_and_unfitted():
    net = _ToyNet(n_layers=3, arch="toy", solver="mt1d")
    assert "unfitted" in repr(net)
    net.fit(
        np.zeros((3, 2), dtype=np.float32),
        np.ones((3, 3), dtype=np.float32),
    )
    assert "fitted" in repr(net)
    assert "toy" in repr(net)


def test_base_emnet_coerce_input_ndarray_passthrough():
    net = _ToyNet(n_layers=2)
    X = np.zeros((3, 2), dtype=np.float64)
    out = net._coerce_input(X)
    assert out.dtype == np.float32


def test_base_emnet_coerce_input_forward_response():
    from pycsamt.forward.em1d import ForwardResponse

    net = _ToyNet(n_layers=2)
    resp = ForwardResponse(
        method="MT1D",
        freqs=np.array([1.0, 10.0, 100.0]),
        rho_a=np.array([10.0, 20.0, 30.0]),
        phase=np.array([40.0, 45.0, 50.0]),
    )
    out = net._coerce_input(resp)
    assert out.shape[0] == 1
    assert out.dtype == np.float32


def test_base_emnet_coerce_input_invalid_type_raises():
    net = _ToyNet(n_layers=2)
    with pytest.raises(TypeError, match="Cannot coerce"):
        net._coerce_input(object())


def test_base_emnet_resolve_device_no_crash():
    net = _ToyNet(n_layers=2, device=None)
    dev = net._resolve_device()
    assert isinstance(dev, str)


# ─────────────────────────────────────────────────────────────────────────────
# BaseEMProcessor — minimal concrete subclass
# ─────────────────────────────────────────────────────────────────────────────


class _ToyProcessor(BaseEMProcessor):
    def __init__(self, scale: float = 1.0):
        self.scale = float(scale)
        self._is_fitted = False

    def fit(self, X, **kwargs):
        self._is_fitted = True
        return self

    def transform(self, X):
        return np.asarray(X) * self.scale

    def _get_params(self):
        return {"scale": self.scale}

    def _get_weights(self):
        return {"scale_w": np.array([self.scale])}

    def _load_weights(self, weights):
        self.scale = float(weights["scale_w"][0])


def test_base_emprocessor_fit_transform():
    proc = _ToyProcessor(scale=2.0)
    X = np.array([1.0, 2.0, 3.0])
    out = proc.fit_transform(X)
    np.testing.assert_allclose(out, [2.0, 4.0, 6.0])
    assert proc._is_fitted


def test_base_emprocessor_save_load_round_trip(tmp_path):
    proc = _ToyProcessor(scale=3.5)
    proc.fit(np.zeros(3))
    path = tmp_path / "toy_proc.npz"
    proc.save(path)

    loaded = _ToyProcessor.load(path)
    assert loaded.scale == pytest.approx(3.5)
    np.testing.assert_allclose(loaded.transform([1.0, 2.0]), [3.5, 7.0])


def test_base_emprocessor_load_pretrained_not_implemented():
    with pytest.raises(NotImplementedError):
        _ToyProcessor.load_pretrained("qc-v1")


def test_base_emprocessor_repr():
    proc = _ToyProcessor()
    assert repr(proc) == "_ToyProcessor()"


def test_base_emprocessor_default_get_weights_and_load_weights_are_noops():
    class _BareProcessor(BaseEMProcessor):
        def fit(self, X, **kwargs):
            return self

        def transform(self, X):
            return X

    proc = _BareProcessor()
    assert proc._get_weights() == {}
    assert proc._get_params() == {}
    proc._load_weights({"anything": np.array([1])})  # must not raise


# ─────────────────────────────────────────────────────────────────────────────
# EMCheckpoint — extra branch not covered by test_inversion.py
# ─────────────────────────────────────────────────────────────────────────────


def test_emcheckpoint_defaults_when_optional_args_omitted(tmp_path):
    ckpt = EMCheckpoint(params={"a": 1})
    assert ckpt.weights == {}
    assert ckpt.history == {}
    assert ckpt.meta == {}

    path = tmp_path / "ckpt_defaults.npz"
    ckpt.save(path)
    loaded = EMCheckpoint.load(path)
    assert loaded.params == {"a": 1}
    assert loaded.weights == {}
