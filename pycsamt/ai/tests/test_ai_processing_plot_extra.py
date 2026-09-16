# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Coverage-gap tests for :mod:`pycsamt.ai.processing.plot`.

Complements ``test_ai_processing_coverage.py`` (qc/denoise/anomaly/classify
extraction logic + ``plot_qc_feature_heatmap``/``plot_qc_summary``) and
``test_ai_processing_api_contracts.py`` (a handful of QC-plot axis-reuse and
input-validation checks). This file targets every plotting function in
``plot.py`` that was previously untested end to end: the shared
``plot_training_history`` helper and the per-estimator families
(denoise / anomaly / dimensionality / imputer / uncertainty / distortion /
time-series-denoise), plus the flexible-input and branch-level gaps left in
the QC family (``_normalise_qc_input``, ``_resolve_profile_colors``,
``_render_bar_chart``, ``_station_freq_grid``).
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pycsamt.ai.plot import StationTickConfig
from pycsamt.ai.processing.plot import (
    plot_anomaly_score_distribution,
    plot_anomaly_scores,
    plot_anomaly_summary,
    plot_denoise_noise_reduction,
    plot_denoise_spectra,
    plot_denoise_summary,
    plot_dimensionality_map,
    plot_dimensionality_summary,
    plot_distortion_feature_space,
    plot_distortion_map,
    plot_distortion_summary,
    plot_imputer_gaps,
    plot_imputer_reconstruction,
    plot_imputer_summary,
    plot_imputer_validation,
    plot_predicted_strike_rose,
    plot_qc_heatmap,
    plot_qc_score_distribution,
    plot_qc_score_spread,
    plot_qc_scores,
    plot_qc_summary,
    plot_training_history,
    plot_ts_denoise_mmf_split,
    plot_ts_denoise_segments,
    plot_ts_denoise_summary,
    plot_uncertainty_map,
    plot_uncertainty_summary,
    plot_uncertainty_validation,
)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# Shared synthetic-data builders
# ─────────────────────────────────────────────────────────────────────────────


def _score_df(
    n_st=6,
    n_fr=4,
    with_profile=False,
    with_station=True,
    with_nan_score=False,
    seed=1,
):
    """A QC-style long table with configurable optional columns, so every
    branch of ``_normalise_qc_input``'s DataFrame handling can be driven
    independently."""
    rng = np.random.default_rng(seed)
    stations = [f"S{i:02d}" for i in range(n_st)]
    freqs = np.array([1000.0, 100.0, 10.0, 1.0][:n_fr])
    rows = []
    for si, st in enumerate(stations):
        for fr in freqs:
            row = {"freq": fr, "score": float(rng.uniform(0, 1))}
            if with_station:
                row["station"] = st
            if with_profile:
                row["profile"] = "L18" if si < n_st // 2 else "L22"
            rows.append(row)
    df = pd.DataFrame(rows)
    if with_nan_score:
        df.loc[df.index[0], "score"] = np.nan
    return df


def _history(n_epochs=12, with_val=True, seed=0):
    rng = np.random.default_rng(seed)
    train = list(np.abs(rng.standard_normal(n_epochs)).cumsum()[::-1] + 0.1)
    hist = {"train_loss": train}
    if with_val:
        hist["val_loss"] = list(np.asarray(train) * 1.15 + 0.02)
    return hist


class _FakeEstimator:
    def __init__(self, history_):
        self.history_ = history_


def _denoise_arrays(n_sites=5, n_comp=4, n_freqs=10, seed=2, rougher=False):
    rng = np.random.default_rng(seed)
    freq = np.logspace(3, -1, n_freqs)
    X_raw = rng.standard_normal((n_sites, n_comp, n_freqs)).astype(float)
    # Denoised = smoothed raw -> strictly lower roughness (positive % reduction)
    kernel = np.array([0.25, 0.5, 0.25])
    X_denoised = np.apply_along_axis(
        lambda v: np.convolve(v, kernel, mode="same"), -1, X_raw
    )
    if rougher:
        # Deliberately make "denoised" rougher than raw for one site so
        # plot_denoise_noise_reduction exercises its negative-color branch.
        X_denoised[0] = X_raw[0] + rng.standard_normal((n_comp, n_freqs)) * 5.0
    return freq, X_raw, X_denoised


def _dim_table(n_st=6, n_fr=4, seed=3, with_nan=False):
    rng = np.random.default_rng(seed)
    stations = [f"S{i:02d}" for i in range(n_st)]
    freqs = np.array([1000.0, 100.0, 10.0, 1.0][:n_fr])
    rows = []
    for st in stations:
        for fr in freqs:
            dim = int(rng.integers(0, 3))
            strike = float(rng.uniform(-90, 90)) if dim == 1 else np.nan
            rows.append(
                {"station": st, "freq": fr, "dim": dim, "strike": strike}
            )
    df = pd.DataFrame(rows)
    if with_nan:
        df.loc[df.index[0], "dim"] = np.nan
    return df


def _gap_table(n_st=5, n_fr=4, seed=4):
    rng = np.random.default_rng(seed)
    stations = [f"S{i:02d}" for i in range(n_st)]
    freqs = np.array([1000.0, 100.0, 10.0, 1.0][:n_fr])
    rows = []
    for st in stations:
        for fr in freqs:
            rows.append(
                {
                    "station": st,
                    "freq": fr,
                    "missing": int(rng.uniform() < 0.25),
                }
            )
    return pd.DataFrame(rows)


def _uncertainty_table(n_st=6, n_fr=4, seed=5):
    rng = np.random.default_rng(seed)
    stations = [f"S{i:02d}" for i in range(n_st)]
    freqs = np.array([1000.0, 100.0, 10.0, 1.0][:n_fr])
    rows = []
    for st in stations:
        for fr in freqs:
            raw = float(rng.uniform(0.01, 0.5))
            rows.append(
                {
                    "station": st,
                    "freq": fr,
                    "z_err_frac": raw,
                    "z_err_frac_calibrated": raw * float(rng.uniform(0.6, 1.0)),
                }
            )
    return pd.DataFrame(rows)


def _distortion_table(n_st=12, seed=6, regimes=(0, 1, 2)):
    rng = np.random.default_rng(seed)
    stations = [f"S{i:02d}" for i in range(n_st)]
    regime = rng.choice(list(regimes), size=n_st)
    df = pd.DataFrame(
        {
            "station": stations,
            "regime": regime,
            "confidence": rng.uniform(0.5, 1.0, n_st),
            "delta_log10_rho": rng.normal(0, 0.3, n_st),
            "twist_deg": rng.uniform(-30, 30, n_st),
        }
    )
    return df


def _ts_denoise_data(n=400, seed=7, with_noisy=True):
    rng = np.random.default_rng(seed)
    t = np.arange(n) * 0.05
    low = np.sin(2 * np.pi * 0.2 * t)
    high = 0.1 * rng.standard_normal(n)
    x = low + high
    denoised = low.copy()
    if with_noisy:
        diagnostics = pd.DataFrame(
            {
                "channel": ["hx"] * 2,
                "win_start": [t[10], t[200]],
                "win_stop": [t[40], t[260]],
                "score": [0.2, 0.1],
                "label": ["noisy", "noisy"],
            }
        )
    else:
        diagnostics = pd.DataFrame(
            columns=["channel", "win_start", "win_stop", "score", "label"]
        )
    return t, x, low, high, denoised, diagnostics


# ─────────────────────────────────────────────────────────────────────────────
# plot_training_history
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_training_history_train_and_val_marks_best():
    hist = _history()
    ax = plot_training_history(hist, title="Training")
    lines = ax.get_lines()
    # 2 loss curves + 1 best-epoch vline
    assert len(lines) == 3
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert "Train" in labels and "Validation" in labels
    assert ax.get_title() == "Training"


def test_plot_training_history_log_y_and_no_mark_best():
    hist = _history()
    ax = plot_training_history(hist, log_y=True, mark_best=False)
    assert ax.get_yscale() == "log"
    assert len(ax.get_lines()) == 2  # no best-epoch vline


def test_plot_training_history_train_only_no_legend():
    hist = _history(with_val=False)
    ax = plot_training_history(hist, show_legend=False)
    assert len(ax.get_lines()) == 1
    assert ax.get_legend() is None


def test_plot_training_history_accepts_fitted_estimator():
    est = _FakeEstimator(_history())
    ax = plot_training_history(est)
    assert len(ax.get_lines()) == 3


def test_plot_training_history_rejects_bad_type():
    with pytest.raises(TypeError, match="history_"):
        plot_training_history(42)


def test_plot_training_history_rejects_empty_history():
    with pytest.raises(ValueError, match="Empty training history"):
        plot_training_history({"train_loss": [], "val_loss": []})
    with pytest.raises(ValueError, match="Empty training history"):
        plot_training_history(_FakeEstimator({}))


# ─────────────────────────────────────────────────────────────────────────────
# EMDenoiser plots
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_denoise_spectra_default_and_explicit_sites():
    freq, X_raw, X_den = _denoise_arrays(n_sites=6, n_comp=4, n_freqs=12)
    fig = plot_denoise_spectra(freq, X_raw, X_den, suptitle="Spectra")
    assert fig.axes  # n_comp * n_show axes
    assert len(fig.axes) == 4 * 4  # default n_show=4
    assert fig.legends  # shared "Input"/"Denoised" legend

    fig2 = plot_denoise_spectra(
        freq, X_raw, X_den, sites=[0, 2, 5], sharey="none",
        component_labels=["a", "b", "c", "d"],
    )
    assert len(fig2.axes) == 4 * 3
    assert fig2.axes[0].get_ylabel() == "a"


def test_plot_denoise_spectra_8_and_arbitrary_component_counts():
    freq, X_raw, X_den = _denoise_arrays(n_sites=2, n_comp=8, n_freqs=8)
    fig = plot_denoise_spectra(freq, X_raw, X_den, n_show=2)
    assert len(fig.axes) == 8 * 2

    freq3, X_raw3, X_den3 = _denoise_arrays(n_sites=2, n_comp=3, n_freqs=8)
    fig3 = plot_denoise_spectra(freq3, X_raw3, X_den3, n_show=1)
    assert fig3.axes[0].get_ylabel() == "ch0"


def test_plot_denoise_noise_reduction_mean_and_median_and_negative():
    freq, X_raw, X_den = _denoise_arrays(rougher=True)
    ax = plot_denoise_noise_reduction(X_raw, X_den, aggregate="mean")
    # first site was made deliberately rougher after denoising -> negative bar
    heights = [p.get_height() for p in ax.patches]
    assert any(h < 0 for h in heights)

    ax2 = plot_denoise_noise_reduction(
        X_raw, X_den, aggregate="median", show_grid=False,
        station_tick_config=StationTickConfig(every=2),
    )
    assert ax2.get_ylabel() == "Roughness reduction (%)"


def test_plot_denoise_summary_with_and_without_history():
    freq, X_raw, X_den = _denoise_arrays(n_sites=4)
    fig = plot_denoise_summary(
        freq, X_raw, X_den, history=_history(), suptitle="Denoiser"
    )
    assert len(fig.axes) >= 3

    fig2 = plot_denoise_summary(freq, X_raw, X_den, history=None)
    assert len(fig2.axes) >= 2

    # a history object whose loss lists are both empty must not add a panel
    fig3 = plot_denoise_summary(
        freq, X_raw, X_den, history=_FakeEstimator({"train_loss": [], "val_loss": []})
    )
    assert len(fig3.axes) == len(fig2.axes)


# ─────────────────────────────────────────────────────────────────────────────
# AnomalyDetector plots
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_anomaly_scores_ndarray_dict_and_dataframe():
    ax = plot_anomaly_scores(np.array([0.1, 0.9, 0.4, 0.2]))
    assert ax.get_legend() is None  # single unnamed profile -> no legend

    scores = {"L18": np.array([0.1, 0.2, 0.9]), "L22": np.array([0.3, 0.8])}
    ax2 = plot_anomaly_scores(
        scores, threshold=0.5, ylim=(0, 1.2), title="Anomaly scores",
        profile_colors=["#111111", "#222222"],
    )
    assert ax2.get_title() == "Anomaly scores"
    assert ax2.get_legend() is not None

    df = pd.DataFrame(
        {"station": ["A", "B", "C"], "score": [0.1, 0.9, 0.3]}
    )
    ax3 = plot_anomaly_scores(df, threshold_percentile=90.0)
    assert ax3.get_ylim()[0] == 0.0


def test_plot_anomaly_score_distribution_kde_and_no_kde():
    scores = {"L18": np.array([0.1, 0.2, 0.9, 0.4]), "L22": np.array([0.3, 0.8])}
    ax = plot_anomaly_score_distribution(scores, show_kde=True)
    assert ax.get_xlabel() == "Anomaly score"

    ax2 = plot_anomaly_score_distribution(
        scores, show_kde=False, show_legend=False, threshold=0.5,
    )
    assert ax2.get_legend() is None


def test_plot_anomaly_summary_two_panels():
    scores = {"L18": np.array([0.1, 0.2, 0.9, 0.4]), "L22": np.array([0.3, 0.8])}
    fig = plot_anomaly_summary(
        scores, suptitle="Anomaly summary",
        station_tick_config=StationTickConfig(every=1),
    )
    assert len(fig.axes) == 2
    assert fig._suptitle.get_text() == "Anomaly summary"


# ─────────────────────────────────────────────────────────────────────────────
# DimensionalityClassifier plots
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_dimensionality_map_validates_and_renders():
    with pytest.raises(ValueError, match="DataFrame must have columns"):
        plot_dimensionality_map(pd.DataFrame({"station": ["S00"]}))

    df = _dim_table(with_nan=True)
    ax = plot_dimensionality_map(df, title="Dim map")
    assert ax.get_ylabel() == "Period (s)"
    assert ax.get_legend() is not None

    ax2 = plot_dimensionality_map(df, station_markers=False, title="No markers")
    assert ax2.get_xlabel() == "Station"


def test_plot_predicted_strike_rose_empty_full_and_dataframe():
    ax_empty = plot_predicted_strike_rose(np.array([np.nan, np.nan]))
    assert ax_empty.get_title() == "No 2-D strike estimates"

    strike = np.array([10.0, 15.0, 170.0, 20.0, np.nan])
    ax = plot_predicted_strike_rose(strike, title="Rose", fold_180=True)
    assert ax.get_title() == "Rose"
    assert len(ax.patches) == 18  # default bins

    ax2 = plot_predicted_strike_rose(strike, fold_180=False, show_mean=False)
    assert len(ax2.patches) == 18

    df = pd.DataFrame({"strike": strike})
    ax3 = plot_predicted_strike_rose(df)
    assert len(ax3.patches) == 18


def test_plot_dimensionality_summary_full_and_missing_strike_col():
    df = _dim_table()
    fig = plot_dimensionality_summary(df, suptitle="Dim summary")
    assert len(fig.axes) == 3

    fig2 = plot_dimensionality_summary(df, strike_col="does_not_exist")
    assert len(fig2.axes) == 3


# ─────────────────────────────────────────────────────────────────────────────
# EMImputer plots
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_imputer_gaps_validates_and_renders():
    with pytest.raises(TypeError, match="DataFrame"):
        plot_imputer_gaps(np.zeros((2, 2)))

    with pytest.raises(ValueError, match="DataFrame must have columns"):
        plot_imputer_gaps(pd.DataFrame({"station": ["S00"]}))

    gt = _gap_table()
    ax = plot_imputer_gaps(gt, title="Gaps")
    assert ax.get_ylabel() == "Period (s)"

    ax2 = plot_imputer_gaps(gt, station_markers=False)
    assert ax2.get_xlabel() == "Station"


def test_plot_imputer_validation_single_panel_and_faceted():
    rng = np.random.default_rng(8)
    y_true = rng.standard_normal(60)
    y_pred = y_true + 0.1 * rng.standard_normal(60)

    fig = plot_imputer_validation(y_true, y_pred, suptitle="Validation")
    assert len(fig.axes) == 1
    assert "RMSE" in fig.axes[0].texts[0].get_text()

    # 3 components with n_cols=2 -> 4 grid slots, 1 hidden
    cids = np.repeat([0, 1, 2], 20)
    fig2 = plot_imputer_validation(
        y_true, y_pred, component_ids=cids, n_cols=2,
        component_labels=["logZxy"], colors=["#123456", "#abcdef", "#00ff00"],
    )
    visible = [a for a in fig2.axes if a.get_visible()]
    hidden = [a for a in fig2.axes if not a.get_visible()]
    assert len(visible) == 3
    assert len(hidden) == 1
    titles = [a.get_title() for a in visible]
    assert "logZxy" in titles
    assert "component 1" in titles  # label list shorter than n components


def test_plot_imputer_validation_handles_too_few_points():
    fig = plot_imputer_validation(np.array([1.0]), np.array([1.0]))
    ax = fig.axes[0]
    # show_metrics is skipped for <=1 point -- no RMSE annotation
    assert not any("RMSE" in t.get_text() for t in ax.texts)


def test_plot_imputer_reconstruction_infers_and_accepts_mask():
    rng = np.random.default_rng(9)
    n_sites, n_comp, n_f = 5, 4, 10
    freq = np.logspace(3, -1, n_f)
    X_raw = rng.standard_normal((n_sites, n_comp, n_f))
    X_raw[0, 0, 2] = np.nan
    X_raw[3, 1, 5] = np.nan
    X_filled = np.nan_to_num(X_raw, nan=0.0)

    fig = plot_imputer_reconstruction(freq, X_raw, X_filled, n_show=3)
    assert len(fig.axes) == n_comp * 3
    assert fig.legends

    mask = np.isnan(X_raw)
    fig2 = plot_imputer_reconstruction(
        freq, X_raw, X_filled, missing_mask=mask, sites=[0, 3], sharey="none",
    )
    assert len(fig2.axes) == n_comp * 2


def test_plot_imputer_summary_with_and_without_history():
    rng = np.random.default_rng(10)
    n_sites, n_comp, n_f = 4, 4, 8
    freq = np.logspace(3, -1, n_f)
    X_raw = rng.standard_normal((n_sites, n_comp, n_f))
    X_raw[1, 0, 3] = np.nan
    X_filled = np.nan_to_num(X_raw, nan=0.0)
    gt = _gap_table(n_st=4, n_fr=4)

    fig = plot_imputer_summary(
        gt, freq, X_raw, X_filled, history=_history(), n_show=2,
        suptitle="Imputer",
    )
    assert len(fig.axes) >= 3

    fig2 = plot_imputer_summary(gt, freq, X_raw, X_filled, history=None)
    assert len(fig2.axes) >= 2


# ─────────────────────────────────────────────────────────────────────────────
# UncertaintyCalibrator plots
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_uncertainty_map_validates_and_renders():
    with pytest.raises(TypeError, match="DataFrame"):
        plot_uncertainty_map(np.zeros((2, 2)))

    with pytest.raises(ValueError, match="DataFrame must have columns"):
        plot_uncertainty_map(pd.DataFrame({"station": ["S00"]}))

    table = _uncertainty_table()
    ax = plot_uncertainty_map(table, title="Uncertainty")
    assert ax.get_ylabel() == "Period (s)"

    ax2 = plot_uncertainty_map(
        table, value_col="z_err_frac", vmin=0.01, vmax=0.5,
        station_markers=False,
    )
    assert ax2.get_xlabel() == "Station"


def test_plot_uncertainty_map_all_zero_values_uses_default_scale():
    table = _uncertainty_table()
    table["z_err_frac_calibrated"] = 0.0
    ax = plot_uncertainty_map(table)
    assert ax is not None


def test_plot_uncertainty_validation_scatter():
    rng = np.random.default_rng(11)
    y_true = rng.uniform(0.01, 0.5, 40)
    y_pred = y_true * rng.uniform(0.8, 1.2, 40)
    ax = plot_uncertainty_validation(y_true, y_pred, title="Held-out")
    assert ax.get_title() == "Held-out"


def test_plot_uncertainty_summary_with_and_without_history():
    table = _uncertainty_table()
    rng = np.random.default_rng(12)
    y_true = rng.uniform(0.01, 0.5, 40)
    y_pred = y_true * rng.uniform(0.8, 1.2, 40)

    fig = plot_uncertainty_summary(
        table, y_true, y_pred, history=_history(), suptitle="Uncertainty"
    )
    assert len(fig.axes) >= 4

    fig2 = plot_uncertainty_summary(table, y_true, y_pred, history=None)
    assert len(fig2.axes) >= 3


# ─────────────────────────────────────────────────────────────────────────────
# DistortionTypeClassifier plots
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_distortion_map_validates_and_renders_both_order_branches():
    with pytest.raises(TypeError, match="DataFrame"):
        plot_distortion_map(np.zeros((2, 2)))

    with pytest.raises(ValueError, match="DataFrame must have columns"):
        plot_distortion_map(pd.DataFrame({"station": ["S00"]}))

    table = _distortion_table()
    ax = plot_distortion_map(table, title="Distortion")
    assert ax.get_ylim() == (0.0, 1.05)

    order = list(table["station"])[::-1]
    ax2 = plot_distortion_map(table, station_labels=order, show_grid=False)
    assert ax2.get_ylabel() == "Confidence"


def test_plot_distortion_feature_space_thresholds_and_missing_class():
    table = _distortion_table(regimes=(0, 2))  # class 1 entirely absent
    ax = plot_distortion_feature_space(table, title="Feature space")
    assert ax.get_xlabel() == "delta_log10_rho"

    with pytest.raises(TypeError, match="DataFrame"):
        plot_distortion_feature_space(np.zeros((2, 2)))
    with pytest.raises(ValueError, match="DataFrame must have columns"):
        plot_distortion_feature_space(pd.DataFrame({"regime": [0]}))

    ax2 = plot_distortion_feature_space(
        table, shift_th=None, twist_th=None, show_legend=False,
        xlabel="dRho", ylabel="Twist",
    )
    assert ax2.get_xlabel() == "dRho"
    assert ax2.get_legend() is None


def test_plot_distortion_summary_with_and_without_history():
    table = _distortion_table()
    fig = plot_distortion_summary(table, history=_history(), suptitle="Distortion")
    assert len(fig.axes) >= 4

    fig2 = plot_distortion_summary(table, history=None)
    assert len(fig2.axes) >= 3


# ─────────────────────────────────────────────────────────────────────────────
# Time-series denoising
# ─────────────────────────────────────────────────────────────────────────────


def test_plot_ts_denoise_mmf_split():
    t, x, low, high, _denoised, _diag = _ts_denoise_data()
    fig = plot_ts_denoise_mmf_split(t, x, low, high, suptitle="MMF split")
    assert len(fig.axes) == 2
    assert fig.axes[0].get_title().startswith("(a)")


def test_plot_ts_denoise_segments_with_and_without_noisy_windows():
    t, x, low, high, _denoised, diag = _ts_denoise_data(with_noisy=True)
    ax = plot_ts_denoise_segments(t, high, diag, title="Residual")
    assert len(ax.patches) >= 2  # shaded noisy spans
    assert ax.get_legend() is not None

    _, _, _, high2, _, empty_diag = _ts_denoise_data(with_noisy=False)
    ax2 = plot_ts_denoise_segments(t, high2, empty_diag)
    assert len(ax2.patches) == 0

    ax3 = plot_ts_denoise_segments(t, high, None)
    assert len(ax3.patches) == 0


def test_plot_ts_denoise_summary_default_and_with_low_high():
    t, x, low, high, denoised, diag = _ts_denoise_data(with_noisy=True)
    fig = plot_ts_denoise_summary(t, x, denoised, diag, suptitle="TS summary")
    assert len(fig.axes) == 3

    fig2 = plot_ts_denoise_summary(t, x, denoised, diag, low=low, high=high)
    assert len(fig2.axes) == 3


# ─────────────────────────────────────────────────────────────────────────────
# QC-family flexible-input & branch gaps
# (`_normalise_qc_input`, `_resolve_profile_colors`, `_render_bar_chart`,
# `_station_freq_grid`, `plot_qc_heatmap`, `plot_qc_score_distribution`,
# `plot_qc_score_spread`)
# ─────────────────────────────────────────────────────────────────────────────


def test_normalise_qc_input_dataframe_profile_without_station():
    df = _score_df(with_profile=True, with_station=False)
    ax = plot_qc_scores(df)
    assert ax.get_legend() is not None  # two profiles -> legend


def test_normalise_qc_input_dataframe_station_without_profile():
    df = _score_df(with_profile=False, with_station=True)
    ax = plot_qc_scores(df)
    assert ax is not None


def test_normalise_qc_input_dataframe_bare_score_column():
    df = _score_df(with_profile=False, with_station=False)
    ax = plot_qc_scores(df)
    assert ax.get_legend() is None  # single "_" profile -> no legend


def test_normalise_qc_input_ndarray_direct():
    ax = plot_qc_scores(np.array([0.2, 0.6, 0.9, 0.1]))
    assert ax is not None


def test_resolve_profile_colors_list_and_dict_forms():
    scores = {"L18": np.array([0.2, 0.5]), "L22": np.array([0.7, 0.9])}
    ax = plot_qc_scores(scores, profile_colors=["#101010", "#202020"])
    patches = ax.get_legend().get_texts()
    assert len(patches) >= 2

    ax2 = plot_qc_scores(
        scores, profile_colors={"L18": "#ff0000", "L22": "#00ff00"}
    )
    assert ax2.get_legend() is not None


def test_render_bar_chart_threshold_lw_zero_skips_line_and_scatter():
    scores = {"L18": np.array([0.2, 0.5, 0.7]), "L22": np.array([0.4, 0.9])}
    ax = plot_qc_scores(
        scores, threshold_lw=0.0, title="No threshold line",
        ylim=(0.0, 1.2),
    )
    assert ax.get_title() == "No threshold line"
    # legend still built (profile patches) but no threshold Line2D entry
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert not any("threshold" in lbl.lower() for lbl in labels)


def test_render_bar_chart_show_scatter_with_profile_column():
    df = _score_df(n_st=4, n_fr=3, with_profile=True, with_station=True)
    ax = plot_qc_scores(df, show_scatter=True)
    assert ax is not None


def test_make_tick_config_override_used_verbatim():
    scores = np.array([0.1, 0.5, 0.9, 0.3, 0.6])
    cfg = StationTickConfig(every=2, rotation=30, fontsize=9)
    ax = plot_qc_scores(scores, station_tick_config=cfg)
    assert ax is not None


def test_qc_score_distribution_drops_empty_profile_and_kde_import_error(
    monkeypatch,
):
    scores = {
        "empty": np.array([np.nan, np.nan]),
        "L18": np.array([0.2, 0.5, 0.8, 0.9, 0.1]),
    }
    ax = plot_qc_score_distribution(scores, show_kde=True)
    assert ax is not None

    monkeypatch.setitem(sys.modules, "scipy.stats", None)
    ax2 = plot_qc_score_distribution(
        {"L18": np.array([0.2, 0.5, 0.8, 0.9, 0.1, 0.3])}, show_kde=True,
    )
    assert ax2 is not None


def test_qc_score_spread_violin_default_and_drops_empty_profile():
    scores = {
        "empty": np.array([np.nan]),
        "L18": np.array([0.2, 0.5, 0.8, 0.9, 0.1]),
        "L22": np.array([0.3, 0.4, 0.6]),
    }
    ax = plot_qc_score_spread(scores)  # default kind="violin"
    assert len(ax.get_xticklabels()) == 2  # empty profile dropped


def test_qc_heatmap_default_render_and_station_markers_false():
    df = _score_df(n_st=5, n_fr=4)
    ax = plot_qc_heatmap(df, title="Heatmap")
    assert ax.get_ylabel() == "Period (s)"

    ax2 = plot_qc_heatmap(df, station_markers=False, show_threshold_contour=True)
    assert ax2.get_xlabel() == "Station"


def test_station_freq_grid_skips_nan_and_unmatched_rows():
    df = _score_df(n_st=4, n_fr=3, with_nan_score=True)
    extra = pd.DataFrame(
        [{"station": "OUTSIDE", "freq": 1000.0, "score": 0.5}]
    )
    df = pd.concat([df, extra], ignore_index=True)
    ax = plot_qc_heatmap(df)
    assert ax is not None


def test_qc_summary_violin_default_kind():
    scores = {"L18": np.array([0.2, 0.5, 0.8]), "L22": np.array([0.3, 0.9])}
    fig = plot_qc_summary(scores, suptitle="QC")
    assert len(fig.axes) >= 3


# ─────────────────────────────────────────────────────────────────────────────
# Remaining narrow branch gaps (explicit station_labels / period_up=False /
# figsize / log_freq toggles, degenerate single-value grids, kde ImportError
# in the anomaly family, and the estimator-with-history_ branch in the
# summary figures, which the dict-history calls above never exercise).
# ─────────────────────────────────────────────────────────────────────────────


def test_qc_heatmap_explicit_station_labels_and_no_contour():
    df = _score_df(n_st=4, n_fr=3)
    order = sorted(df["station"].unique())
    ax = plot_qc_heatmap(df, station_labels=order, show_threshold_contour=False)
    assert ax is not None


def test_qc_feature_heatmap_remaining_branches():
    from pycsamt.ai.processing.plot import plot_qc_feature_heatmap

    df = _score_df(n_st=5, n_fr=3)
    df["snr"] = np.linspace(2, 20, len(df))
    order = sorted(df["station"].unique())

    fig = plot_qc_feature_heatmap(
        df, features=["snr"], station_labels=order, figsize=(7.0, 3.0),
        period_up=False, station_markers=False, title="Feature heatmap",
    )
    assert len(fig.axes) >= 1

    # single-frequency grid (n_f == 1 branch) + an entirely-NaN feature
    # column (fin.size == 0 -> default vmin/vmax branch) + a row outside
    # the station order (skip branch inside the fill loop).
    single_freq = df[df["freq"] == df["freq"].iloc[0]].copy()
    single_freq["snr"] = np.nan
    extra = single_freq.iloc[[0]].copy()
    extra["station"] = "OUTSIDE"
    single_freq = pd.concat([single_freq, extra], ignore_index=True)
    fig2 = plot_qc_feature_heatmap(single_freq, features=["snr"])
    assert len(fig2.axes) >= 1


def test_denoise_spectra_log_freq_false():
    freq, X_raw, X_den = _denoise_arrays(n_sites=2, n_comp=2, n_freqs=6)
    fig = plot_denoise_spectra(freq, X_raw, X_den, log_freq=False, n_show=2)
    assert fig.axes[0].get_xscale() == "linear"


def test_resolve_anomaly_threshold_all_nan_pooled_scores():
    scores = {"L18": np.array([np.nan, np.nan])}
    ax = plot_anomaly_scores(scores)
    assert ax is not None


def test_anomaly_score_distribution_kde_import_error(monkeypatch):
    monkeypatch.setitem(sys.modules, "scipy.stats", None)
    scores = np.array([0.1, 0.2, 0.3, 0.9, 0.8, 0.4])
    ax = plot_anomaly_score_distribution(scores, show_kde=True)
    assert ax is not None


def test_dimensionality_map_period_up_false_and_title_no_markers():
    df = _dim_table()
    ax = plot_dimensionality_map(df, period_up=False)
    assert ax is not None

    ax2 = plot_dimensionality_map(
        df, station_markers=False, title="No markers title"
    )
    assert ax2.get_title() == "No markers title"


def test_imputer_gaps_explicit_station_labels_and_period_up_false():
    gt = _gap_table()
    order = sorted(gt["station"].unique())
    ax = plot_imputer_gaps(gt, station_labels=order, period_up=False)
    assert ax is not None


def test_validation_scatter_panel_show_metrics_false():
    rng = np.random.default_rng(20)
    y_true = rng.standard_normal(30)
    y_pred = y_true + 0.05 * rng.standard_normal(30)
    fig = plot_imputer_validation(y_true, y_pred, show_metrics=False)
    ax = fig.axes[0]
    assert not any("RMSE" in t.get_text() for t in ax.texts)


def test_imputer_validation_explicit_figsize_both_modes():
    rng = np.random.default_rng(21)
    y_true = rng.standard_normal(40)
    y_pred = y_true + 0.05 * rng.standard_normal(40)

    fig = plot_imputer_validation(y_true, y_pred, figsize=(4.0, 4.0))
    assert fig.get_size_inches().tolist() == [4.0, 4.0]

    cids = np.repeat([0, 1], 20)
    fig2 = plot_imputer_validation(
        y_true, y_pred, component_ids=cids, figsize=(6.0, 3.0)
    )
    assert fig2.get_size_inches().tolist() == [6.0, 3.0]


def test_imputer_reconstruction_log_freq_false_and_no_missing_values():
    rng = np.random.default_rng(22)
    n_sites, n_comp, n_f = 3, 2, 6
    freq = np.logspace(2, 0, n_f)
    X_raw = rng.standard_normal((n_sites, n_comp, n_f))  # no NaNs at all
    X_filled = X_raw.copy()

    fig = plot_imputer_reconstruction(freq, X_raw, X_filled, log_freq=False)
    assert fig.axes[0].get_xscale() == "linear"


def test_uncertainty_map_explicit_station_labels_and_period_up_false():
    table = _uncertainty_table()
    order = sorted(table["station"].unique())
    ax = plot_uncertainty_map(table, station_labels=order, period_up=False)
    assert ax is not None


def test_summary_figures_accept_fitted_estimator_history():
    table = _distortion_table()
    est = _FakeEstimator(_history())
    fig = plot_distortion_summary(table, history=est)
    assert len(fig.axes) >= 4

    utable = _uncertainty_table()
    rng = np.random.default_rng(23)
    y_true = rng.uniform(0.01, 0.5, 20)
    y_pred = y_true * rng.uniform(0.8, 1.2, 20)
    fig2 = plot_uncertainty_summary(
        utable, y_true, y_pred, history=_FakeEstimator(_history())
    )
    assert len(fig2.axes) >= 4


def test_ts_denoise_functions_without_suptitle():
    t, x, low, high, denoised, diag = _ts_denoise_data()
    fig = plot_ts_denoise_mmf_split(t, x, low, high)
    assert fig._suptitle is None


def test_dimensionality_map_explicit_station_labels():
    df = _dim_table()
    order = sorted(df["station"].unique())
    ax = plot_dimensionality_map(df, station_labels=order)
    assert ax is not None


def test_imputer_reconstruction_with_suptitle():
    rng = np.random.default_rng(24)
    n_sites, n_comp, n_f = 2, 2, 6
    freq = np.logspace(2, 0, n_f)
    X_raw = rng.standard_normal((n_sites, n_comp, n_f))
    X_raw[0, 0, 1] = np.nan
    X_filled = np.nan_to_num(X_raw, nan=0.0)
    fig = plot_imputer_reconstruction(
        freq, X_raw, X_filled, suptitle="Reconstruction"
    )
    assert fig._suptitle.get_text() == "Reconstruction"


def test_imputer_summary_explicit_missing_mask():
    rng = np.random.default_rng(25)
    n_sites, n_comp, n_f = 3, 2, 6
    freq = np.logspace(2, 0, n_f)
    X_raw = rng.standard_normal((n_sites, n_comp, n_f))
    X_raw[1, 0, 2] = np.nan
    X_filled = np.nan_to_num(X_raw, nan=0.0)
    mask = np.isnan(X_raw)
    gt = _gap_table(n_st=3, n_fr=4)
    fig = plot_imputer_summary(gt, freq, X_raw, X_filled, missing_mask=mask)
    assert len(fig.axes) >= 2
