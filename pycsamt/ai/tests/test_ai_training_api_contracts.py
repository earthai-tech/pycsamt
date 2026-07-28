"""API-contract tests for :mod:`pycsamt.ai.training`.

The older training tests cover common workflows.  These tests focus on
contract details: public exports, dataset filtering/inversion, normalizer edge
cases, deterministic augmentation, metric NaN behavior, and trainer state.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from pycsamt.ai.training import (
    AugmentFreqDrop,
    AugmentMixup,
    AugmentNoise,
    AugmentStaticShift,
    Compose,
    EMDataset,
    EMTrainer,
    Normalizer,
    RandomApply,
    depth_rmse,
    layer_rmse,
    mae,
    masked_mse_loss,
    r2,
    relative_rmse,
    rmse,
    summarise,
)


def _forward_dataset() -> SimpleNamespace:
    X = np.arange(30, dtype=np.float32).reshape(6, 5)
    y = np.array(
        [
            [1.0, 2.0, 3.0, 10.0, 100.0, np.nan, np.nan],
            [2.0, 3.0, 4.0, 20.0, 200.0, np.nan, np.nan],
            [3.0, 4.0, 5.0, 30.0, 300.0, np.nan, np.nan],
            [4.0, 5.0, 6.0, 40.0, 400.0, 700.0, 900.0],
            [5.0, 6.0, 7.0, 50.0, 500.0, 800.0, 1000.0],
            [6.0, 7.0, 8.0, 60.0, 600.0, 900.0, 1100.0],
        ],
        dtype=np.float32,
    )
    meta = {"n_layers": np.array([3, 3, 3, 4, 4, 4])}
    return SimpleNamespace(X=X, y=y, meta=meta)


def test_ai_training_package_exports_public_api():
    assert Normalizer.__name__ == "Normalizer"
    assert EMDataset.__name__ == "EMDataset"
    assert EMTrainer.__name__ == "EMTrainer"
    assert callable(rmse)
    assert callable(mae)
    assert callable(r2)
    assert callable(relative_rmse)
    assert callable(depth_rmse)
    assert callable(layer_rmse)
    assert callable(masked_mse_loss)
    assert callable(summarise)
    assert AugmentNoise.__name__ == "AugmentNoise"
    assert AugmentStaticShift.__name__ == "AugmentStaticShift"
    assert AugmentFreqDrop.__name__ == "AugmentFreqDrop"
    assert AugmentMixup.__name__ == "AugmentMixup"
    assert Compose.__name__ == "Compose"
    assert RandomApply.__name__ == "RandomApply"


def test_normalizer_unfitted_errors_and_zero_variance_roundtrip():
    norm = Normalizer(eps=1e-6)
    X = np.ones((4, 3), dtype=np.float32)

    with pytest.raises(RuntimeError, match=r"fit\(\) before transform"):
        norm.transform(X)
    with pytest.raises(RuntimeError, match=r"fit\(\) before inverse_transform"):
        norm.inverse_transform(X)

    Xn = norm.fit_transform(X)
    assert np.all(np.isfinite(Xn))
    assert np.allclose(Xn, 0.0)
    np.testing.assert_allclose(norm.inverse_transform(Xn), X)
    assert np.all(norm.std > 0.0)


def test_normalizer_dict_roundtrip_preserves_eps_independent_statistics():
    X = np.array([[1.0, np.nan], [3.0, 5.0], [5.0, 9.0]])
    norm = Normalizer(eps=1e-4).fit(X)
    restored = Normalizer.from_dict(norm.to_dict())

    np.testing.assert_allclose(restored.mean, norm.mean)
    np.testing.assert_allclose(restored.std, norm.std)
    np.testing.assert_allclose(restored.transform(X), norm.transform(X))


def test_emdataset_filters_layers_logs_thickness_and_inverse_roundtrip():
    fds = _forward_dataset()
    ds = EMDataset(fds, n_layers=3, log_thickness=True)

    assert len(ds) == 3
    assert ds.n_features == 5
    assert ds.n_params == 5
    assert ds._n_layers == 3
    assert "EMDataset" in repr(ds)

    y_raw = ds.inverse_y(ds.y)
    np.testing.assert_allclose(y_raw[:, :3], fds.y[:3, :3], atol=1e-5)
    np.testing.assert_allclose(y_raw[:, 3:], fds.y[:3, 3:5], atol=1e-4)
    np.testing.assert_allclose(ds.inverse_x(ds.X), fds.X[:3], atol=1e-5)


def test_emdataset_split_shares_normalizers_and_disables_val_augmentation():
    ds = EMDataset(_forward_dataset(), n_layers=4, augment_noise=0.2)
    train, val = ds.split(val_frac=1 / 3, seed=0)

    assert len(train) + len(val) == len(ds)
    assert train.x_norm is val.x_norm is ds.x_norm
    assert train.y_norm is val.y_norm is ds.y_norm
    assert train.augment_noise == pytest.approx(0.2)
    assert val.augment_noise == 0.0


def test_emdataset_rejects_missing_layer_count():
    with pytest.raises(ValueError, match="No samples with n_layers=9"):
        EMDataset(_forward_dataset(), n_layers=9)


def test_augmenters_are_deterministic_with_explicit_rng_and_preserve_y():
    X = np.ones((5, 8), dtype=np.float32)
    y = np.arange(10, dtype=np.float32).reshape(5, 2)
    rng1 = np.random.default_rng(42)
    rng2 = np.random.default_rng(42)
    aug = Compose(
        [
            AugmentStaticShift(0.1, n_amp_features=4, per_sample=True),
            AugmentNoise(sigma=0.02, clip=1.0),
            AugmentFreqDrop(drop_rate=0.25, contiguous=True, fill_value=-9.0),
        ],
        seed=0,
    )

    Xa1, ya1 = aug.augmenters[0](X, y, rng=rng1)
    Xa2, ya2 = aug.augmenters[0](X, y, rng=rng2)

    np.testing.assert_allclose(Xa1, Xa2)
    assert ya1 is y
    assert ya2 is y

    Xc, yc = aug(X, y)
    assert Xc.shape == X.shape
    assert yc is y
    assert np.any(Xc == -9.0)


def test_random_apply_and_mixup_edge_contracts():
    X = np.arange(12, dtype=np.float32).reshape(3, 4)
    y = np.arange(6, dtype=np.float32).reshape(3, 2)

    never = RandomApply(AugmentNoise(sigma=10.0), p=0.0, seed=0)
    X_never, y_never = never(X, y)
    np.testing.assert_array_equal(X_never, X)
    assert y_never is y
    assert "RandomApply" in repr(never)

    single = X[:1]
    single_y = y[:1]
    X_single, y_single = AugmentMixup(alpha=0.5)(single, single_y)
    assert X_single is single
    assert y_single is single_y


def test_metrics_ignore_nan_and_report_all_nan_as_nan():
    y_true = np.array([[1.0, 2.0, np.nan], [2.0, 4.0, 6.0]])
    y_pred = np.array([[1.5, 1.0, 3.0], [2.5, 5.0, 7.0]])

    assert rmse(y_true, y_pred) > 0.0
    assert mae(y_true, y_pred) > 0.0
    assert np.isfinite(r2(y_true, y_pred))
    assert relative_rmse(y_true, y_pred) > 0.0
    assert depth_rmse(y_true, y_pred, n_layers=2) > 0.0
    assert layer_rmse(y_true, y_pred).shape == (3,)

    summary = summarise(y_true, y_pred, n_layers=2)
    assert set(summary) == {
        "rmse",
        "mae",
        "r2",
        "relative_rmse",
        "depth_rmse",
    }

    all_nan = np.full((2, 2), np.nan)
    assert np.isnan(rmse(all_nan, all_nan))
    assert np.isnan(mae(all_nan, all_nan))
    assert np.isnan(r2(all_nan, all_nan))
    assert np.isnan(relative_rmse(all_nan, all_nan))
    assert np.isnan(depth_rmse(all_nan, all_nan, n_layers=2))


def test_emtrainer_initial_state_and_repr():
    model = object()
    trainer = EMTrainer(
        model,
        lr=2e-3,
        weight_decay=0.0,
        patience=3,
        min_delta=1e-4,
        batch_size=16,
        device="cpu",
        grad_clip=1.0,
        verbose=False,
    )

    assert trainer.model is model
    assert trainer.history == {
        "train_loss": [],
        "val_loss": [],
        "lr": [],
        "epoch_time": [],
    }
    assert trainer.best_epoch == 0
    assert trainer.best_val_loss == float("inf")
    assert "EMTrainer" in repr(trainer)
    assert "batch_size=16" in repr(trainer)


def test_masked_mse_loss_uses_only_finite_targets():
    torch = pytest.importorskip("torch")

    pred = torch.tensor([[1.0, 2.0], [3.0, 4.0]], requires_grad=True)
    target = torch.tensor([[1.0, float("nan")], [1.0, 6.0]])

    loss = masked_mse_loss(pred, target)
    # three finite targets: (1-1)^2, (3-1)^2 and (4-6)^2 -> mean over 3
    assert float(loss.detach()) == pytest.approx(
        (0.0 + (3.0 - 1.0) ** 2 + (4.0 - 6.0) ** 2) / 3
    )

    zero = masked_mse_loss(pred, torch.full_like(target, float("nan")))
    assert float(zero.detach()) == 0.0


def test_emtrainer_tiny_torch_fit_and_weight_roundtrip():
    torch = pytest.importorskip("torch")

    class TinyDataset:
        def __init__(self, X, y):
            self.X = X.astype(np.float32)
            self.y = y.astype(np.float32)

        def __len__(self):
            return len(self.X)

        def __getitem__(self, idx):
            return torch.from_numpy(self.X[idx]), torch.from_numpy(self.y[idx])

    rng = np.random.default_rng(0)
    X = rng.normal(size=(12, 3)).astype(np.float32)
    true_w = np.array([[1.0], [-2.0], [0.5]], dtype=np.float32)
    y = X @ true_w
    train = TinyDataset(X[:9], y[:9])
    val = TinyDataset(X[9:], y[9:])
    model = torch.nn.Linear(3, 1)

    trainer = EMTrainer(
        model,
        lr=1e-2,
        batch_size=4,
        patience=5,
        device="cpu",
        verbose=False,
    )
    trainer.fit(train, val, epochs=2)

    assert len(trainer.history["train_loss"]) >= 1
    assert len(trainer.history["val_loss"]) >= 1
    assert np.isfinite(trainer.best_val_loss)

    weights = trainer.get_weights()
    clone = torch.nn.Linear(3, 1)
    clone_trainer = EMTrainer(clone, verbose=False)
    clone_trainer.load_weights(weights)

    for name, value in clone.state_dict().items():
        np.testing.assert_allclose(value.detach().numpy(), weights[name])
