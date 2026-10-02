# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Additional coverage for pycsamt.ai.processing.denoise.

test_ai_processing_coverage.py already exercises EMDenoiser's torch
fit/transform/save/load round trip, verbose training, and the scipy
Gaussian-fallback path. This file fills the remaining gap: input
validation errors, calling transform()/apply() before fit(), the full
sites-in/sites-out ``apply()`` path (including the frequency-grid
interpolation branch in ``prepare_z_features``/``_reconstruct_z_block``
and ``n_components=8``), the ``ImportError`` branches when emtools is
unavailable, the "no valid Z data" branch, and ``__repr__``/``history_``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from pycsamt.ai.processing.denoise import EMDenoiser, prepare_z_features

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


def _denoise_X(n=24, n_comp=4, n_freqs=16, seed=5):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n, n_comp, n_freqs)).astype(np.float32)


# ---------------------------------------------------------------------------
# fit() validation
# ---------------------------------------------------------------------------


def test_fit_rejects_non_3d_input():
    den = EMDenoiser()
    with pytest.raises(ValueError, match="must be 3-D"):
        den.fit(np.zeros((5, 4)), epochs=1)


def test_fit_rejects_mismatched_n_freqs_on_second_call():
    den = EMDenoiser(n_freqs=10)
    with pytest.raises(ValueError, match="Expected n_freqs=10"):
        den.fit(_denoise_X(n=4, n_freqs=8), epochs=1)


# ---------------------------------------------------------------------------
# Calling before fit()
# ---------------------------------------------------------------------------


def test_transform_before_fit_raises():
    den = EMDenoiser()
    with pytest.raises(RuntimeError, match="Call fit"):
        den.transform(_denoise_X(n=2))


def test_apply_before_fit_raises():
    den = EMDenoiser()
    with pytest.raises(RuntimeError, match="Call fit"):
        den.apply(object())


# ---------------------------------------------------------------------------
# __repr__ / history_
# ---------------------------------------------------------------------------


def test_repr_unfitted_and_history_empty():
    den = EMDenoiser()
    assert repr(den) == "EMDenoiser(n_freqs=None, n_components=4, unfitted)"
    assert den.history_ == {}


def test_repr_fitted():
    den = EMDenoiser(channels=(4, 8, 4))
    den.fit(_denoise_X(n=8, n_freqs=10), epochs=1, verbose=False)
    assert "fitted" in repr(den)
    assert den.history_["train_loss"]


# ---------------------------------------------------------------------------
# prepare_z_features
# ---------------------------------------------------------------------------


def test_prepare_z_features_import_error(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def _fake_import(name, *a, **kw):
        if name == "pycsamt.emtools._core":
            raise ImportError("blocked")
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", _fake_import)
    with pytest.raises(ImportError, match="emtools is required"):
        prepare_z_features([1, 2, 3])


def test_prepare_z_features_no_valid_z_raises(monkeypatch):
    import pycsamt.ai.processing.denoise as denoise_mod

    class _FakeCore:
        @staticmethod
        def ensure_sites(sites, **kw):
            return [object()]

        @staticmethod
        def _iter_items(sites):
            return iter(sites)

        @staticmethod
        def _get_z_block(ed, with_errors=False):
            return (None, None, None)

    monkeypatch.setitem(sys.modules, "pycsamt.emtools._core", _FakeCore)
    with pytest.raises(ValueError, match="No valid Z data"):
        prepare_z_features([object()])


def test_interp_channel_to_grid_recovers_sampled_points():
    from pycsamt.ai.processing.denoise import _interp_channel_to_grid

    # Descending frequency grid, as real EDI/MT data is conventionally
    # stored -- exercises the internal ascending-sort-before-interp step.
    src_freq = np.array([100.0, 10.0, 1.0, 0.1])
    values = np.array([4.0, 3.0, 2.0, 1.0])
    out = _interp_channel_to_grid(values, src_freq, src_freq)
    np.testing.assert_allclose(out, values)


def test_interp_channel_to_grid_onto_a_different_grid():
    from pycsamt.ai.processing.denoise import _interp_channel_to_grid

    src_freq = np.array([100.0, 10.0, 1.0])
    values = np.array([3.0, 2.0, 1.0])
    dst_freq = np.array([50.0, 5.0])
    out = _interp_channel_to_grid(values, src_freq, dst_freq)
    assert out.shape == dst_freq.shape
    assert np.all(np.isfinite(out))


def _pack_feat(z, n_components):
    """Mirror prepare_z_features's own packing, without ensure_sites."""

    def _amp(zc):
        return np.log10(np.maximum(np.abs(zc), 1e-24))

    def _phase(zc):
        return np.degrees(np.angle(zc))

    comps = [_amp(z[:, 0, 1]), _phase(z[:, 0, 1]), _amp(z[:, 1, 0]), _phase(z[:, 1, 0])]
    if n_components == 8:
        comps += [_amp(z[:, 0, 0]), _phase(z[:, 0, 0]), _amp(z[:, 1, 1]), _phase(z[:, 1, 1])]
    return np.stack(comps, axis=0)


def test_reconstruct_z_block_round_trip_4_components():
    from pycsamt.ai.processing.denoise import _reconstruct_z_block

    freq = np.logspace(2, -1, 8)
    z = np.zeros((8, 2, 2), dtype=complex)
    z[:, 0, 1] = 10.0 + 5.0j
    z[:, 1, 0] = -8.0 - 2.0j
    z[:, 0, 0] = 1.0 + 1.0j  # untouched diagonal

    feat = _pack_feat(z, 4)
    z2 = _reconstruct_z_block(z, freq, freq, feat, 4, log_amp=True)
    np.testing.assert_allclose(z2[:, 0, 1], z[:, 0, 1], rtol=1e-4)
    np.testing.assert_allclose(z2[:, 1, 0], z[:, 1, 0], rtol=1e-4)
    np.testing.assert_array_equal(z2[:, 0, 0], z[:, 0, 0])


def test_reconstruct_z_block_round_trip_8_components():
    from pycsamt.ai.processing.denoise import _reconstruct_z_block

    freq = np.logspace(2, -1, 6)
    z = np.zeros((6, 2, 2), dtype=complex)
    z[:, 0, 1] = 10.0 + 5.0j
    z[:, 1, 0] = -8.0 - 2.0j
    z[:, 0, 0] = 1.0 + 0.5j
    z[:, 1, 1] = 0.9 - 0.3j

    feat = _pack_feat(z, 8)
    z2 = _reconstruct_z_block(z, freq, freq, feat, 8, log_amp=True)
    np.testing.assert_allclose(z2[:, 0, 0], z[:, 0, 0], rtol=1e-4)
    np.testing.assert_allclose(z2[:, 1, 1], z[:, 1, 1], rtol=1e-4)
    np.testing.assert_allclose(z2[:, 0, 1], z[:, 0, 1], rtol=1e-4)
    np.testing.assert_allclose(z2[:, 1, 0], z[:, 1, 0], rtol=1e-4)


# ---------------------------------------------------------------------------
# apply() — full sites-in / sites-out round trip
# ---------------------------------------------------------------------------


def test_apply_round_trip_4_components(sites):
    from pycsamt.emtools._core import _get_z_block, _iter_items, ensure_sites

    den = EMDenoiser(channels=(4, 8, 4))
    X = prepare_z_features(sites, n_components=4)
    den.fit(X, epochs=2, verbose=False)

    corrected = den.apply(sites, inplace=False)
    assert corrected is not None

    orig_first = next(_iter_items(ensure_sites(sites)))
    corr_first = next(_iter_items(ensure_sites(corrected)))
    _, z_orig, _ = _get_z_block(orig_first, with_errors=False)[:3]
    _, z_corr, _ = _get_z_block(corr_first, with_errors=False)[:3]
    assert z_corr.shape == z_orig.shape
    assert np.all(np.isfinite(z_corr[:, 0, 1]))
    assert np.all(np.isfinite(z_corr[:, 1, 0]))


def test_apply_round_trip_8_components_inplace(sites):
    den = EMDenoiser(n_components=8, channels=(4, 8, 4))
    X = prepare_z_features(sites, n_components=8)
    den.fit(X, epochs=2, verbose=False)

    corrected = den.apply(sites, inplace=True)
    assert corrected is not None


def test_apply_import_error(monkeypatch):
    import builtins

    den = EMDenoiser()
    den._is_fitted = True  # bypass the "call fit" guard to reach the import

    real_import = builtins.__import__

    def _fake_import(name, *a, **kw):
        if name == "pycsamt.emtools._core":
            raise ImportError("blocked")
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", _fake_import)
    with pytest.raises(ImportError, match="emtools is required for"):
        den.apply(object())
