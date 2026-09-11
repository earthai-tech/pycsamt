# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for :mod:`pycsamt.ai.processing.tsdenoise` (MMF-SVM-K-SVD)."""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.ai.processing import (
    KSVDDenoiser,
    SignalQualityClassifier,
    TimeSeriesDenoiser,
    compute_entropy_features,
    generate_synthetic_library,
    mmf_split,
    plot_ts_denoise_mmf_split,
    plot_ts_denoise_segments,
    plot_ts_denoise_summary,
    snr_db,
    time_domain_ncc,
)
from pycsamt.ai.processing.tsdenoise import (
    _synth_clean,
    _synth_pulse,
    _synth_square,
    approx_entropy,
    box_dimension,
    fuzzy_entropy,
    ksvd_dictionary,
    omp,
    omp_batch,
    sample_entropy,
)
from pycsamt.ts import TSData

# --------------------------------------------------------------- metrics


def test_snr_and_ncc_identical_and_degenerate():
    rng = np.random.default_rng(0)
    x = rng.normal(size=200)
    assert snr_db(x, x) == float("inf")
    assert time_domain_ncc(x, x) == pytest.approx(1.0, abs=1e-9)
    assert time_domain_ncc(x, -x) == pytest.approx(-1.0, abs=1e-9)
    assert np.isnan(time_domain_ncc(np.zeros(10), np.zeros(10)))


# --------------------------------------------------------------- MMF


def test_mmf_split_reconstructs_input():
    rng = np.random.default_rng(1)
    x = rng.normal(size=400)
    low, high = mmf_split(x, size=15)
    assert np.allclose(low + high, x, atol=1e-8)


def test_mmf_split_requires_size_or_dt():
    x = np.zeros(50)
    with pytest.raises(ValueError):
        mmf_split(x)


def test_mmf_split_nan_preserved():
    x = np.arange(100, dtype=float)
    x[10:13] = np.nan
    low, high = mmf_split(x, size=9)
    assert np.all(np.isnan(low[10:13]))
    assert np.all(np.isnan(high[10:13]))
    assert not np.any(np.isnan(low[:10]))


def test_mmf_split_size_must_exceed_pulse_width_to_reject_it():
    """A pulse wider than the structuring element leaks into `low`;
    widening `size` past the pulse width rejects it into `high`."""
    rng = np.random.default_rng(2)
    n = 400
    x = _synth_clean(n, rng, amp=1.0)
    x[150:190] += 40.0  # 40-sample-wide plateau

    low_narrow, _ = mmf_split(x, size=9)
    low_wide, _ = mmf_split(x, size=81)

    assert low_narrow[165] > 5.0  # plateau leaks through a narrow SE
    assert abs(low_wide[165]) < 5.0  # a wide-enough SE rejects it


# --------------------------------------------------------------- entropy


def test_entropy_features_nan_on_constant_or_short_input():
    const = np.ones(20)
    assert np.isnan(sample_entropy(const))
    assert np.isnan(fuzzy_entropy(const))
    assert np.isnan(approx_entropy(const))
    too_short = np.array([1.0, 2.0])
    assert np.isnan(sample_entropy(too_short))
    assert np.isnan(box_dimension(too_short))


def test_entropy_features_separate_clean_from_noisy():
    """High-quality (random) segments should score higher SE/FE/AE
    than a stereotyped, low-complexity pulse -- the separation the
    SVM classifier relies on (Fig. 4 of Gui et al., 2024)."""
    rng = np.random.default_rng(3)
    clean = _synth_clean(150, rng, amp=1.0)
    pulse = _synth_pulse(150, rng, amp=8.0)

    fc = compute_entropy_features(clean[None, :])[0]
    fp = compute_entropy_features(pulse[None, :])[0]

    assert np.all(np.isfinite(fc))
    assert np.all(np.isfinite(fp))
    assert fc[0] > fp[0]  # sample entropy
    assert fc[1] > fp[1]  # fuzzy entropy


def test_compute_entropy_features_batch_shape():
    rng = np.random.default_rng(4)
    segs = rng.normal(size=(12, 100))
    F = compute_entropy_features(segs)
    assert F.shape == (12, 4)


# ------------------------------------------------- synthetic library


def test_generate_synthetic_library_shapes_and_labels():
    segs, y = generate_synthetic_library(120, n_per_class=15, seed=0)
    assert segs.shape == (30, 120)
    assert set(np.unique(y)) == {0, 1}
    assert (y == 1).sum() == 15
    assert (y == 0).sum() == 15


def test_generate_synthetic_library_deterministic_with_seed():
    a, ya = generate_synthetic_library(80, n_per_class=10, seed=7)
    b, yb = generate_synthetic_library(80, n_per_class=10, seed=7)
    assert np.allclose(a, b)
    assert np.array_equal(ya, yb)


# ------------------------------------------- SignalQualityClassifier


def test_signal_quality_classifier_holdout_accuracy():
    segs, y = generate_synthetic_library(150, n_per_class=150, seed=11)
    rng = np.random.default_rng(12)
    idx = rng.permutation(len(segs))
    n_train = int(0.7 * len(segs))
    tr, te = idx[:n_train], idx[n_train:]

    clf = SignalQualityClassifier(random_state=0).fit(segs[tr], y[tr])
    pred = clf.predict(segs[te])
    acc = (pred == y[te]).mean()
    assert acc > 0.85  # paper reports a clean linear separation (Fig. 6)


def test_signal_quality_classifier_from_synthetic_factory():
    clf = SignalQualityClassifier.from_synthetic(
        win_len=150, n_per_class=40, random_state=0,
    )
    assert clf._is_fitted
    scores = clf.transform(
        generate_synthetic_library(150, n_per_class=5, seed=1)[0]
    )
    assert scores.shape == (10,)
    assert np.all((scores >= 0.0) & (scores <= 1.0))


def test_signal_quality_classifier_save_load_roundtrip(tmp_path):
    clf = SignalQualityClassifier.from_synthetic(
        win_len=120, n_per_class=30, random_state=0,
    )
    segs, _ = generate_synthetic_library(120, n_per_class=5, seed=2)
    before = clf.transform(segs)

    path = tmp_path / "sqc.npz"
    clf.save(path)
    loaded = SignalQualityClassifier.load(path)
    after = loaded.transform(segs)

    assert np.allclose(before, after)


def test_signal_quality_classifier_requires_fit_before_transform():
    clf = SignalQualityClassifier()
    with pytest.raises(RuntimeError):
        clf.transform(np.zeros((2, 50)))


# --------------------------------------------------------------- OMP / K-SVD


def test_omp_recovers_sparse_code_from_known_dictionary():
    rng = np.random.default_rng(5)
    L, K = 20, 8
    D = rng.normal(size=(L, K))
    D /= np.linalg.norm(D, axis=0)
    true_code = np.zeros(K)
    true_code[[2, 5]] = [3.0, -1.5]
    y = D @ true_code

    code = omp(D, y, sparsity=2)
    assert np.allclose(code, true_code, atol=1e-6)


def test_omp_batch_shape():
    rng = np.random.default_rng(6)
    D = rng.normal(size=(16, 10))
    D /= np.linalg.norm(D, axis=0)
    Y = rng.normal(size=(16, 5))
    X = omp_batch(D, Y, sparsity=3)
    assert X.shape == (10, 5)


def test_ksvd_dictionary_reduces_reconstruction_error():
    rng = np.random.default_rng(7)
    L, N = 16, 60
    Y = rng.normal(size=(L, N)) * 0.1
    # inject a repeated, low-complexity pattern most columns share
    pattern = np.sin(np.linspace(0, 3.14, L))
    Y[:, ::2] += pattern[:, None] * 3.0

    D, X = ksvd_dictionary(Y, n_atoms=6, sparsity=2, n_iter=8, seed=0)
    assert D.shape == (L, 6)
    assert np.allclose(np.linalg.norm(D, axis=0), 1.0, atol=1e-6)

    err_learned = np.linalg.norm(Y - D @ X)
    # a random, unlearned dictionary of the same size should do worse
    D_rand = rng.normal(size=(L, 6))
    D_rand /= np.linalg.norm(D_rand, axis=0)
    X_rand = omp_batch(D_rand, Y, sparsity=2)
    err_random = np.linalg.norm(Y - D_rand @ X_rand)

    assert err_learned < err_random


def test_ksvd_denoiser_reduces_transient_noise():
    rng = np.random.default_rng(8)
    n = 200
    background = _synth_clean(n, rng, amp=1.0)
    noise = _synth_square(n, rng, amp=10.0)
    noisy = background + noise

    den = KSVDDenoiser(
        patch_len=16, n_atoms=16, sparsity=2, n_iter=10, random_state=0,
    ).fit(noisy)
    denoised = den.transform(noisy)

    assert snr_db(background, denoised) > snr_db(background, noisy)


def test_ksvd_denoiser_save_load_roundtrip(tmp_path):
    rng = np.random.default_rng(9)
    x = _synth_clean(150, rng) + _synth_pulse(150, rng, amp=6.0)
    den = KSVDDenoiser(patch_len=16, n_atoms=12, sparsity=2, random_state=0)
    den.fit(x)
    before = den.transform(x)

    path = tmp_path / "ksvd.npz"
    den.save(path)
    loaded = KSVDDenoiser.load(path)
    after = loaded.transform(x)

    assert np.allclose(before, after)


def test_ksvd_denoiser_requires_fit_before_noise_profile():
    den = KSVDDenoiser()
    with pytest.raises(RuntimeError):
        den.noise_profile(np.zeros(50))


# ------------------------------------------------- TimeSeriesDenoiser


def _make_noisy_ts(seed: int = 0, n: int = 3000, dt: float = 1 / 15.0):
    rng = np.random.default_rng(seed)
    ex = _synth_clean(n, rng, amp=1.0)
    ex[1000:1010] += 25.0
    ex[2000:2016] -= 30.0
    ey = _synth_clean(n, rng, amp=1.0)
    return TSData(data={"EX": ex, "EY": ey}, dt=dt, station="synthtest")


def test_time_series_denoiser_apply_roundtrip_shapes():
    ts = _make_noisy_ts()
    den = TimeSeriesDenoiser(mmf_size=61, random_state=0)
    out = den.apply(ts)

    assert set(out.channels()) == set(ts.channels())
    assert out.get("EX").shape == ts.get("EX").shape
    assert out.get("EY").shape == ts.get("EY").shape

    diag = den.diagnostics_
    assert isinstance(diag, pd.DataFrame)
    assert set(diag["channel"].unique()) == {"EX", "EY"}
    assert set(diag["label"].unique()) <= {"clean", "noisy"}


def test_time_series_denoiser_flags_injected_pulses_as_noisy():
    ts = _make_noisy_ts()
    den = TimeSeriesDenoiser(mmf_size=61, win_seconds=10.0, random_state=0)
    den.apply(ts)
    diag = den.diagnostics_
    ex_diag = diag[diag["channel"] == "EX"]

    t_pulse1 = 1000 * (1 / 15.0)
    t_pulse2 = 2000 * (1 / 15.0)
    hit1 = ex_diag[
        (ex_diag["win_start"] <= t_pulse1) & (t_pulse1 < ex_diag["win_stop"])
    ]
    hit2 = ex_diag[
        (ex_diag["win_start"] <= t_pulse2) & (t_pulse2 < ex_diag["win_stop"])
    ]
    assert (hit1["label"] == "noisy").all()
    assert (hit2["label"] == "noisy").all()


def test_time_series_denoiser_channels_subset_and_not_inplace():
    ts = _make_noisy_ts()
    original_ex = ts.get("EX").copy()

    den = TimeSeriesDenoiser(mmf_size=61, random_state=0)
    out = den.apply(ts, channels=["EX"], inplace=False)

    assert np.allclose(ts.get("EX"), original_ex)  # untouched
    assert np.allclose(out.get("EY"), ts.get("EY"))  # not requested, copied
    assert not np.allclose(out.get("EX"), original_ex)  # denoised


def test_time_series_denoiser_inplace_mutates():
    ts = _make_noisy_ts()
    original_ex = ts.get("EX").copy()
    den = TimeSeriesDenoiser(mmf_size=61, random_state=0)
    out = den.apply(ts, channels=["EX"], inplace=True)

    assert out is ts
    assert not np.allclose(ts.get("EX"), original_ex)


def test_time_series_denoiser_handles_nan_gaps():
    ts = _make_noisy_ts()
    x = ts.get("EX").copy()
    x[500:505] = np.nan
    ts.add_channel("EX", x)

    den = TimeSeriesDenoiser(mmf_size=61, random_state=0)
    out = den.apply(ts, channels=["EX"])
    assert np.all(np.isnan(out.get("EX")[500:505]))
    assert np.all(np.isfinite(out.get("EX")[:500]))


def test_time_series_denoiser_requires_dt():
    ts = TSData(data={"EX": np.zeros(100)}, station="nodts")
    den = TimeSeriesDenoiser(random_state=0)
    with pytest.raises(ValueError):
        den.apply(ts)


def test_time_series_denoiser_reuses_supplied_classifier():
    clf = SignalQualityClassifier.from_synthetic(
        win_len=150, n_per_class=30, random_state=0,
    )
    den = TimeSeriesDenoiser(quality_clf=clf, mmf_size=61)
    assert den._is_fitted
    ts = _make_noisy_ts()
    out = den.apply(ts, channels=["EX"])
    assert den._clf is clf
    assert out.get("EX").shape == ts.get("EX").shape


def test_time_series_denoiser_save_load_roundtrip(tmp_path):
    ts = _make_noisy_ts()
    den = TimeSeriesDenoiser(mmf_size=61, random_state=0)
    den.fit(dt=ts.dt, n_per_class=30)
    out_before = den.apply(ts, channels=["EX"])

    path = tmp_path / "tsden.npz"
    den.save(path)
    loaded = TimeSeriesDenoiser.load(path)
    out_after = loaded.apply(ts, channels=["EX"])

    assert np.allclose(out_before.get("EX"), out_after.get("EX"))


# --------------------------------------------------------------- plots


def test_plot_ts_denoise_mmf_split_returns_figure():
    rng = np.random.default_rng(10)
    n = 300
    t = np.arange(n) * (1 / 15.0)
    x = _synth_clean(n, rng)
    low, high = mmf_split(x, size=31)
    fig = plot_ts_denoise_mmf_split(t, x, low, high)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_plot_ts_denoise_segments_reuses_supplied_axes():
    rng = np.random.default_rng(11)
    n = 300
    t = np.arange(n) * (1 / 15.0)
    high = _synth_clean(n, rng)
    diag = pd.DataFrame(
        {
            "channel": ["EX"] * 3,
            "win_start": [0.0, 10.0, 20.0],
            "win_stop": [10.0, 20.0, 30.0],
            "score": [0.9, 0.1, 0.95],
            "label": ["clean", "noisy", "clean"],
        }
    )
    fig, ax = plt.subplots()
    out_ax = plot_ts_denoise_segments(t, high, diag, ax=ax)
    assert out_ax is ax
    plt.close(fig)


def test_plot_ts_denoise_summary_builds_three_panels():
    ts = _make_noisy_ts(n=1500)
    den = TimeSeriesDenoiser(mmf_size=61, random_state=0)
    out = den.apply(ts, channels=["EX"])
    diag = den.diagnostics_

    t = ts.time("EX")
    x = ts.get("EX")
    low, high = mmf_split(x, size=61)
    fig = plot_ts_denoise_summary(
        t, x, out.get("EX"), diag, low=low, high=high,
    )
    assert isinstance(fig, plt.Figure)
    assert len(fig.axes) == 3
    plt.close(fig)
