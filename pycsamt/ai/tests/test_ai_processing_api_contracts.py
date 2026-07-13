"""API-contract tests for :mod:`pycsamt.ai.processing`.

The existing processing tests cover many estimator workflows.  This module
adds focused public-API checks that stay deterministic when deep-learning
backends are not installed, or when they are installed but should not be used
for a fast contract suite.
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.ai.processing import (
    AnomalyDetector,
    DimensionalityClassifier,
    EMDenoiser,
    EMQCScorer,
    plot_qc_heatmap,
    plot_qc_score_distribution,
    plot_qc_score_spread,
    plot_qc_scores,
    prepare_z_features,
)
from pycsamt.ai.processing.classify import (
    _df_to_feature_matrix,
    _rule_labels,
)
from pycsamt.ai.processing.qc import _extract_qc_features


def _processing_features(n: int = 24) -> np.ndarray:
    rng = np.random.default_rng(123)
    X = rng.normal(size=(n, 5)).astype(np.float32)
    X[:, 0] = np.abs(X[:, 0]) * 4.0
    X[:, 1] = np.abs(X[:, 1]) * 0.3
    return X


def _qc_dataframe() -> pd.DataFrame:
    stations = ["S00", "S01", "S02"]
    freqs = np.array([100.0, 10.0, 1.0])
    rows = []
    for si, station in enumerate(stations):
        for fi, freq in enumerate(freqs):
            score = 0.25 + 0.1 * si + 0.08 * fi
            rows.append(
                {
                    "profile": "L18",
                    "station": station,
                    "freq": freq,
                    "snr": 5.0 + si + fi,
                    "swift_skew": 0.05 * fi,
                    "asym": -0.1 + 0.1 * si,
                    "phase_xy": 35.0 + fi,
                    "phase_yx": -140.0 + fi,
                    "score": score,
                    "flag": int(score >= 0.5),
                }
            )
    return pd.DataFrame(rows)


def test_ai_processing_package_exports_public_api():
    assert EMDenoiser.__name__ == "EMDenoiser"
    assert AnomalyDetector.__name__ == "AnomalyDetector"
    assert DimensionalityClassifier.__name__ == "DimensionalityClassifier"
    assert EMQCScorer.__name__ == "EMQCScorer"
    assert callable(prepare_z_features)
    assert callable(plot_qc_scores)
    assert callable(plot_qc_heatmap)
    assert callable(plot_qc_score_distribution)
    assert callable(plot_qc_score_spread)


def test_emd_noiser_uses_numpy_fallback_when_no_backend(monkeypatch):
    import pycsamt.ai.processing.denoise as denoise_mod

    monkeypatch.setattr(denoise_mod, "active_backend", lambda: "none")
    rng = np.random.default_rng(10)
    X = rng.normal(size=(8, 4, 12)).astype(np.float32)

    den = EMDenoiser(n_components=4, channels=(4, 8, 4), dropout=0.0)
    den.fit(X, epochs=1, seed=0, verbose=False)
    out = den.transform(X)

    assert den.n_freqs == 12
    assert den._backend_name == "numpy"
    assert den._use_numpy is True
    assert out.shape == X.shape
    assert out.dtype == np.float32
    assert np.all(np.isfinite(out))
    assert "fitted" in repr(den)
    assert den._get_params()["channels"] == [4, 8, 4]


def test_emd_noiser_validates_input_shape_and_frequency_count(monkeypatch):
    import pycsamt.ai.processing.denoise as denoise_mod

    monkeypatch.setattr(denoise_mod, "active_backend", lambda: "none")
    den = EMDenoiser(n_freqs=8)

    with pytest.raises(RuntimeError, match=r"fit\(\) before transform"):
        den.transform(np.ones((2, 4, 8), dtype=np.float32))

    with pytest.raises(ValueError, match="3-D"):
        den.fit(np.ones((2, 8), dtype=np.float32), verbose=False)

    with pytest.raises(ValueError, match="Expected n_freqs=8"):
        den.fit(np.ones((2, 4, 7), dtype=np.float32), verbose=False)


def test_anomaly_detector_pca_fallback_transform_flags_and_save_load(
    monkeypatch,
    tmp_path,
):
    # torch/tf are patched out below to force the sklearn PCA path
    pytest.importorskip("sklearn")
    import pycsamt.ai.processing.anomaly as anomaly_mod

    monkeypatch.setattr(
        anomaly_mod.AnomalyDetector,
        "_fit_torch",
        lambda self, Xn, **kwargs: (_ for _ in ()).throw(
            ImportError("no torch")
        ),
    )
    monkeypatch.setattr(
        anomaly_mod.AnomalyDetector,
        "_fit_tensorflow",
        lambda self, Xn, **kwargs: (_ for _ in ()).throw(
            ImportError("no tf")
        ),
    )
    X = np.vstack(
        [
            np.zeros((8, 6), dtype=np.float32),
            np.ones((8, 6), dtype=np.float32),
            np.eye(6, dtype=np.float32),
        ]
    )

    det = AnomalyDetector(
        n_features=6,
        latent_dim=3,
        channels=(8,),
        threshold_percentile=80.0,
    )
    det.fit(X, epochs=1, verbose=False)
    scores = det.transform(X)
    flags = det.flag_anomalies(X)

    assert det._backend_name == "pca"
    assert det._use_pca is True
    assert scores.shape == (len(X),)
    assert flags.shape == (len(X),)
    assert flags.dtype == bool
    assert np.all(np.isfinite(scores))
    assert "pca" in repr(det)

    path = tmp_path / "anomaly_detector.npz"
    det.save(path)
    loaded = AnomalyDetector.load(path)
    np.testing.assert_allclose(loaded.transform(X), scores)
    np.testing.assert_array_equal(loaded.flag_anomalies(X), flags)


def test_anomaly_detector_rejects_feature_mismatch(monkeypatch):
    import pycsamt.ai.processing.anomaly as anomaly_mod

    monkeypatch.setattr(
        anomaly_mod.AnomalyDetector,
        "_fit_torch",
        lambda self, Xn, **kwargs: (_ for _ in ()).throw(
            ImportError("no torch")
        ),
    )
    det = AnomalyDetector(n_features=4)

    with pytest.raises(ValueError, match="Expected n_features=4"):
        det.fit(np.ones((5, 3), dtype=np.float32), verbose=False)


def test_dimensionality_classifier_rule_helpers_and_dataframe_coercion():
    beta = np.array([0.5, 0.5, 4.0, 6.0])
    ellipt = np.array([0.1, 0.4, 0.1, 0.8])
    labels = _rule_labels(beta, ellipt, skew_th=3.0, ellipt_th=0.2)
    np.testing.assert_array_equal(labels, np.array([0, 1, 2, 2]))

    df = pd.DataFrame(
        {
            "beta_abs": [1.0, np.nan],
            "ellipt_abs": [0.1, 0.3],
            "logrho_det": [2.0, 2.5],
        }
    )
    X = _df_to_feature_matrix(df)

    assert X.shape == (2, 5)
    assert X[0, 0] == 1.0
    assert X[1, 0] == 0.0
    assert X[0, 3] == 0.0


def test_dimensionality_classifier_rf_fallback_and_from_features_table(
    monkeypatch,
):
    # torch/tf are patched out below to force the sklearn RF path
    pytest.importorskip("sklearn")
    import pycsamt.ai.processing.classify as classify_mod

    monkeypatch.setattr(
        classify_mod.DimensionalityClassifier,
        "_fit_torch",
        lambda self, Xn, y, strike, **kwargs: (_ for _ in ()).throw(
            ImportError("no torch")
        ),
    )
    monkeypatch.setattr(
        classify_mod.DimensionalityClassifier,
        "_fit_tensorflow",
        lambda self, Xn, y, strike, **kwargs: (_ for _ in ()).throw(
            ImportError("no tf")
        ),
    )
    X = _processing_features(48)
    y = _rule_labels(X[:, 0], X[:, 1])
    df = pd.DataFrame(
        X,
        columns=[
            "beta_abs",
            "ellipt_abs",
            "logrho_det",
            "phi_det",
            "tip_amp",
        ],
    )
    df["dim"] = y

    clf = DimensionalityClassifier(hidden=(16,), dropout=0.0, n_classes=3)
    clf.fit(df, verbose=False)
    proba = clf.transform(df)
    pred = clf.predict(df)
    strike = clf.predict_strike(df)

    assert clf._backend_name == "rf"
    assert clf._use_rf is True
    assert proba.shape == (len(df), 3)
    np.testing.assert_allclose(proba.sum(axis=1), 1.0, atol=1e-7)
    assert pred.shape == (len(df),)
    assert np.all(np.isnan(strike))
    assert "rf" in repr(clf)

    trained = DimensionalityClassifier.from_features_table(
        df,
        epochs=1,
        verbose=False,
    )
    assert trained._is_fitted is True
    assert trained.predict(df).shape == (len(df),)


def test_emqc_scorer_rule_only_scores_dataframe_and_feature_matrix(tmp_path):
    feat = np.array(
        [
            [10.0, 0.1, 0.0, 30.0, -130.0],
            [1.0, 0.1, 0.0, 30.0, -130.0],
            [10.0, 0.9, 0.0, 30.0, -130.0],
            [np.nan, np.nan, 0.0, 30.0, -130.0],
        ],
        dtype=float,
    )
    df = pd.DataFrame(
        feat,
        columns=["snr", "swift_skew", "asym", "phase_xy", "phase_yx"],
    )
    scorer = EMQCScorer(
        snr_threshold=3.0,
        skew_threshold=0.3,
        score_threshold=0.5,
        use_ml=False,
    )

    assert scorer.fit(feat) is scorer
    scores = scorer.transform(feat)
    df_scores = scorer.transform(df)

    np.testing.assert_allclose(scores, np.array([1.0, 0.0, 0.0, 1.0]))
    np.testing.assert_allclose(df_scores, scores)
    assert "rule_only" in repr(scorer)

    path = tmp_path / "qc_rule_only.npz"
    scorer.save(path)
    loaded = EMQCScorer.load(path)
    np.testing.assert_allclose(loaded.transform(feat), scores)


def test_qc_feature_extraction_from_impedance_tensor():
    z = np.zeros((3, 2, 2), dtype=complex)
    z[:, 0, 1] = np.array([1.0 + 1.0j, 2.0 + 0.0j, 4.0 + 0.0j])
    z[:, 1, 0] = np.array([-1.0 + 1.0j, -2.0 + 0.0j, -4.0 + 0.0j])
    z[:, 0, 0] = 0.05
    z[:, 1, 1] = 0.05
    ze = np.full_like(z, 0.1, dtype=float)

    features = _extract_qc_features(z, ze)

    assert features.shape == (3, 5)
    assert np.all(features[:, 0] > 0.0)
    assert np.all(features[:, 1] >= 0.0)
    assert np.all(np.isfinite(features))


def test_qc_plot_functions_reuse_supplied_axes():
    df = _qc_dataframe()
    scores = {"L18": np.array([0.2, 0.7, 0.9]), "L22": np.array([0.4, 0.6])}

    fig1, ax1 = plt.subplots()
    out_ax1 = plot_qc_scores(
        scores,
        score_threshold=0.5,
        show_scatter=False,
        show_zone_labels=False,
        ax=ax1,
    )
    assert out_ax1 is ax1
    assert ax1.get_ylabel() == "QC score"

    fig2, ax2 = plt.subplots()
    out_ax2 = plot_qc_heatmap(df, ax=ax2, show_threshold_contour=False)
    assert out_ax2 is ax2
    assert ax2.get_xlabel() == "Station"

    fig3, ax3 = plt.subplots()
    out_ax3 = plot_qc_score_distribution(
        df,
        show_kde=False,
        show_legend=False,
        ax=ax3,
    )
    assert out_ax3 is ax3
    assert ax3.get_xlabel() == "QC score"

    fig4, ax4 = plt.subplots()
    out_ax4 = plot_qc_score_spread(
        scores,
        kind="box",
        show_scatter=False,
        ax=ax4,
    )
    assert out_ax4 is ax4
    assert ax4.get_ylabel() == "QC score"

    for fig in (fig1, fig2, fig3, fig4):
        plt.close(fig)


def test_qc_plot_functions_validate_inputs():
    with pytest.raises(TypeError, match="DataFrame"):
        plot_qc_heatmap(np.ones((2, 2)))

    with pytest.raises(ValueError, match="'station', 'freq', 'score'"):
        plot_qc_heatmap(pd.DataFrame({"station": ["S00"], "score": [1.0]}))

    with pytest.raises(TypeError, match="ndarray, dict, or DataFrame"):
        plot_qc_scores([0.1, 0.2])
