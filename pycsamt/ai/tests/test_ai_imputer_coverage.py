# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Coverage tests for pycsamt.ai.processing.imputer.EMImputer.

Mirrors the fixture/test-writing conventions of
``test_ai_processing_coverage.py`` (EMDenoiser et al.): real PyTorch
training paths with small epoch counts, the numpy fallback via a
monkeypatched ``active_backend``, save/load round trips, and the
sites-based ``apply()`` path using the bundled ``data/3edis`` EDI
fixture.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_EDI_DIR = _PROJECT_ROOT / "data" / "3edis"
_HAS_EDIS = _EDI_DIR.exists() and any(_EDI_DIR.glob("*.edi"))


@pytest.fixture(scope="module")
def sites():
    if not _HAS_EDIS:
        pytest.skip(f"3edis dataset not found: {_EDI_DIR}")
    from pycsamt.agents import MTLoaderAgent

    r = MTLoaderAgent().execute({"path": str(_EDI_DIR)})
    if r.status != "success":
        pytest.skip("Could not load 3edis sites.")
    return r["sites"]


def _imputer_X(n=24, n_comp=4, n_freqs=16, gap_frac=0.1, seed=3):
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, n_comp, n_freqs)).astype(np.float32)
    gaps = rng.random(X.shape) < gap_frac
    X[gaps] = np.nan
    return X, gaps


# ─────────────────────────────────────────────────────────────────────────────
# construction / basic properties
# ─────────────────────────────────────────────────────────────────────────────


def test_no_args_init():
    from pycsamt.ai.processing.imputer import EMImputer

    imp = EMImputer()
    assert imp.n_freqs is None
    assert imp.n_components == 4
    assert imp.channels == (64, 128, 64)
    assert imp._is_fitted is False
    assert imp.history_ == {}


def test_repr_before_and_after_fit():
    from pycsamt.ai.processing.imputer import EMImputer

    imp = EMImputer(channels=(4, 8, 4))
    assert "unfitted" in repr(imp)

    X, _ = _imputer_X(n=12, n_freqs=10)
    imp.fit(X, epochs=1, verbose=False)
    assert "fitted" in repr(imp)
    assert "n_freqs=10" in repr(imp)


def test_transform_before_fit_raises():
    from pycsamt.ai.processing.imputer import EMImputer

    imp = EMImputer()
    with pytest.raises(RuntimeError):
        imp.transform(np.zeros((2, 4, 8), dtype=np.float32))


def test_apply_before_fit_raises():
    from pycsamt.ai.processing.imputer import EMImputer

    imp = EMImputer()
    with pytest.raises(RuntimeError):
        imp.apply(object())


# ─────────────────────────────────────────────────────────────────────────────
# fit() validation branches
# ─────────────────────────────────────────────────────────────────────────────


def test_fit_infers_n_freqs():
    from pycsamt.ai.processing.imputer import EMImputer

    X, _ = _imputer_X(n=12, n_freqs=20)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=1, verbose=False)
    assert imp.n_freqs == 20


def test_fit_n_freqs_mismatch_raises():
    from pycsamt.ai.processing.imputer import EMImputer

    X, _ = _imputer_X(n=8, n_freqs=16)
    imp = EMImputer(n_freqs=32)
    with pytest.raises(ValueError, match="Expected n_freqs"):
        imp.fit(X, epochs=1, verbose=False)


def test_fit_requires_3d_input():
    from pycsamt.ai.processing.imputer import EMImputer

    imp = EMImputer()
    with pytest.raises(ValueError, match="3-D"):
        imp.fit(np.zeros((4, 4)), epochs=1)


def test_fit_all_nan_raises():
    from pycsamt.ai.processing.imputer import EMImputer

    X = np.full((4, 4, 8), np.nan, dtype=np.float32)
    imp = EMImputer()
    with pytest.raises(ValueError, match="nothing to learn"):
        imp.fit(X, epochs=1)


# ─────────────────────────────────────────────────────────────────────────────
# torch training path
# ─────────────────────────────────────────────────────────────────────────────


def test_fit_torch_backend_selected():
    from pycsamt.ai.processing.imputer import EMImputer

    X, _ = _imputer_X(n=20, n_freqs=12)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=2, verbose=False)
    assert imp._backend_name == "torch"
    assert imp._use_numpy is False
    assert set(imp.history_) == {"train_loss", "val_loss"}
    assert len(imp.history_["train_loss"]) == 2
    assert len(imp.history_["val_loss"]) == 2


def test_fit_torch_verbose_prints(capsys):
    from pycsamt.ai.processing.imputer import EMImputer

    X, _ = _imputer_X(n=16, n_freqs=10)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=4, verbose=True)
    out = capsys.readouterr().out
    assert "Epoch" in out


def test_transform_preserves_observed_and_fills_missing():
    from pycsamt.ai.processing.imputer import EMImputer

    X, gaps = _imputer_X(n=24, n_freqs=16, gap_frac=0.15)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=3, verbose=False)

    X_filled = imp.transform(X)
    assert X_filled.shape == X.shape
    assert np.all(np.isfinite(X_filled))
    # observed cells are byte-identical to the input
    np.testing.assert_array_equal(X_filled[~gaps], X[~gaps])
    # missing cells are no longer NaN
    assert not np.any(np.isnan(X_filled[gaps]))


def test_fit_explicit_init_args_used():
    from pycsamt.ai.processing.imputer import EMImputer

    X, _ = _imputer_X(n=16, n_freqs=8)
    imp = EMImputer(
        n_freqs=8,
        n_components=4,
        channels=(4, 8, 4),
        dropout=0.2,
        device="cpu",
    )
    imp.fit(
        X,
        mask_frac=0.3,
        epochs=2,
        batch_size=4,
        lr=1e-3,
        val_frac=0.2,
        seed=1,
        verbose=False,
    )
    assert imp._is_fitted is True


# ─────────────────────────────────────────────────────────────────────────────
# save / load round trip
# ─────────────────────────────────────────────────────────────────────────────


def test_save_load_round_trip_torch(tmp_path):
    from pycsamt.ai.processing.imputer import EMImputer

    X, _ = _imputer_X(n=20, n_freqs=12)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=2, verbose=False)
    out_before = imp.transform(X)

    path = tmp_path / "imputer.npz"
    imp.save(path)
    loaded = EMImputer.load(path)

    assert loaded._backend_name == "torch"
    out_after = loaded.transform(X)
    np.testing.assert_allclose(out_before, out_after, rtol=1e-4, atol=1e-5)


# ─────────────────────────────────────────────────────────────────────────────
# numpy fallback (_freq_axis_fill)
# ─────────────────────────────────────────────────────────────────────────────


def test_get_weights_on_unfitted_instance_is_empty():
    from pycsamt.ai.processing.imputer import EMImputer

    imp = EMImputer()
    assert imp._get_weights() == {}


def test_numpy_fallback_save_load_round_trip(monkeypatch, tmp_path):
    import pycsamt.ai.processing.imputer as imputer_mod

    monkeypatch.setattr(imputer_mod, "active_backend", lambda: "none")

    X, gaps = _imputer_X(n=10, n_freqs=10, gap_frac=0.2)
    imp = imputer_mod.EMImputer()
    imp.fit(X, epochs=1, verbose=False)
    assert imp._use_numpy is True

    path = tmp_path / "imputer_numpy.npz"
    imp.save(path)
    loaded = imputer_mod.EMImputer.load(path)

    assert loaded._backend_name == "numpy"
    assert loaded._use_numpy is True
    out = loaded.transform(X)
    np.testing.assert_array_equal(out[~gaps], X[~gaps])


def test_numpy_fallback_via_no_backend(monkeypatch):
    import pycsamt.ai.processing.imputer as imputer_mod

    monkeypatch.setattr(imputer_mod, "active_backend", lambda: "none")

    X, gaps = _imputer_X(n=10, n_freqs=12, gap_frac=0.2)
    imp = imputer_mod.EMImputer()
    imp.fit(X, epochs=1, verbose=True)
    assert imp._use_numpy is True
    assert imp._backend_name == "numpy"
    assert imp.history_ == {}

    out = imp.transform(X)
    assert out.shape == X.shape
    assert np.all(np.isfinite(out))
    np.testing.assert_array_equal(out[~gaps], X[~gaps])


def test_numpy_fallback_prints_message(monkeypatch, capsys):
    import pycsamt.ai.processing.imputer as imputer_mod

    monkeypatch.setattr(imputer_mod, "active_backend", lambda: "none")

    X, _ = _imputer_X(n=8, n_freqs=10)
    imp = imputer_mod.EMImputer()
    imp.fit(X, epochs=1, verbose=True)
    out = capsys.readouterr().out
    assert "fallback" in out.lower()


def test_freq_axis_fill_full_nan_channel_uses_mean():
    from pycsamt.ai.processing.imputer import _freq_axis_fill

    X = np.full((1, 2, 6), np.nan, dtype=np.float32)
    X[0, 1] = np.array([1.0, np.nan, 3.0, np.nan, 5.0, 6.0])
    x_mean = np.array([[[7.0], [0.0]]], dtype=np.float32)  # (1, 2, 1)

    out = _freq_axis_fill(X, x_mean)
    assert np.all(out[0, 0] == 7.0)
    assert np.all(np.isfinite(out[0, 1]))
    np.testing.assert_allclose(out[0, 1, [0, 2, 4, 5]], [1.0, 3.0, 5.0, 6.0])


def test_freq_axis_fill_partial_nan_row_interpolated():
    from pycsamt.ai.processing.imputer import _freq_axis_fill

    X = np.zeros((1, 1, 5), dtype=np.float32)
    X[0, 0] = [0.0, np.nan, 2.0, np.nan, 4.0]
    x_mean = np.zeros((1, 1, 1), dtype=np.float32)

    out = _freq_axis_fill(X, x_mean)
    np.testing.assert_allclose(out[0, 0], [0.0, 1.0, 2.0, 3.0, 4.0])


# ─────────────────────────────────────────────────────────────────────────────
# apply() — sites-in / sites-out
# ─────────────────────────────────────────────────────────────────────────────


def test_apply_fills_missing_rows_without_inplace(sites):
    from pycsamt.ai.processing.denoise import prepare_z_features
    from pycsamt.ai.processing.imputer import EMImputer
    from pycsamt.emtools._core import _get_z_block, _iter_items, ensure_sites

    X = prepare_z_features(sites, n_components=4)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=2, verbose=False)

    S = ensure_sites(sites)
    ed0 = next(_iter_items(S))
    _, z0, _ = _get_z_block(ed0, with_errors=False)[:3]

    filled = imp.apply(sites, inplace=False)
    assert filled is not None

    # original untouched
    S_after = ensure_sites(sites)
    ed0_after = next(_iter_items(S_after))
    _, z0_after, _ = _get_z_block(ed0_after, with_errors=False)[:3]
    np.testing.assert_array_equal(
        np.isnan(z0), np.isnan(z0_after)
    )


def test_apply_reconstructs_a_genuinely_missing_row(sites):
    import copy

    from pycsamt.ai.processing.denoise import prepare_z_features
    from pycsamt.ai.processing.imputer import EMImputer
    from pycsamt.emtools._core import _get_z_block, _iter_items, ensure_sites

    S = ensure_sites(copy.deepcopy(sites))
    ed0 = next(_iter_items(S))
    Z0, z0, fr0 = _get_z_block(ed0, with_errors=False)[:3]
    z0 = z0.copy()
    hole_idx = len(z0) // 2
    z0[hole_idx] = np.nan
    Z0.z = z0

    X = prepare_z_features(S, n_components=4)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=3, verbose=False)

    filled = imp.apply(S, inplace=False)
    S_filled = ensure_sites(filled)
    ed0_filled = next(_iter_items(S_filled))
    _, z0_filled, _ = _get_z_block(ed0_filled, with_errors=False)[:3]

    # with the default n_components=4, only Zxy/Zyx are modelled --
    # the diagonal at the reconstructed row is left exactly as it was
    # (NaN), per EMImputer.apply's documented "not modelled -> untouched"
    # contract.
    assert not np.isnan(z0_filled[hole_idx, 0, 1])
    assert not np.isnan(z0_filled[hole_idx, 1, 0])
    assert np.isnan(z0_filled[hole_idx, 0, 0])
    assert np.isnan(z0_filled[hole_idx, 1, 1])
    # rows away from the hole are untouched
    other = np.arange(len(z0)) != hole_idx
    np.testing.assert_array_equal(z0_filled[other], z0[other])


def test_apply_inplace_true_returns_same_collection_type(sites):
    from pycsamt.ai.processing.denoise import prepare_z_features
    from pycsamt.ai.processing.imputer import EMImputer

    X = prepare_z_features(sites, n_components=4)
    imp = EMImputer(channels=(4, 8, 4))
    imp.fit(X, epochs=1, verbose=False)

    result = imp.apply(sites, inplace=True)
    assert type(result) is type(sites) or hasattr(result, "__iter__")
