# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Focused coverage-gap tests for pycsamt.ai.processing.

Complements test_processing.py and test_ai_processing_api_contracts.py by
exercising branches those files don't reach: the ML/torch training paths
(this environment has PyTorch but not TensorFlow installed), sites-based
feature extraction (via the bundled ``data/3edis`` EDI files), verbose
logging branches, and save/load round trips.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
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


# ─────────────────────────────────────────────────────────────────────────────
# qc.py
# ─────────────────────────────────────────────────────────────────────────────


def test_extract_qc_features_without_errors_uses_smoothness_snr():
    from pycsamt.ai.processing.qc import _extract_qc_features

    rng = np.random.default_rng(0)
    n = 12
    zxy = (1.0 + 0.1 * rng.standard_normal(n)) + 1j * (
        1.0 + 0.1 * rng.standard_normal(n)
    )
    zyx = -zxy + 0.05 * rng.standard_normal(n)
    z = np.zeros((n, 2, 2), dtype=complex)
    z[:, 0, 1] = zxy
    z[:, 1, 0] = zyx
    z[:, 0, 0] = 0.02
    z[:, 1, 1] = 0.02

    F = _extract_qc_features(z, None)

    assert F.shape == (n, 5)
    assert np.all(np.isfinite(F[:, 0]))
    assert np.all(F[:, 0] > 0)


def test_sites_to_feature_df_and_prepare_z_features(sites):
    from pycsamt.ai.processing.denoise import prepare_z_features
    from pycsamt.ai.processing.qc import _sites_to_feature_df

    df = _sites_to_feature_df(sites)
    assert set(
        [
            "station",
            "freq",
            "snr",
            "swift_skew",
            "asym",
            "phase_xy",
            "phase_yx",
        ]
    ).issubset(df.columns)
    assert len(df) > 0

    X4 = prepare_z_features(sites, n_components=4)
    assert X4.ndim == 3 and X4.shape[1] == 4
    assert X4.dtype == np.float32

    X8 = prepare_z_features(sites, n_components=8)
    assert X8.shape[1] == 8


def test_emqc_scorer_fit_transform_with_ml():
    from pycsamt.ai.processing.qc import EMQCScorer

    rng = np.random.default_rng(1)
    feat = np.column_stack(
        [
            rng.uniform(2.0, 20.0, 60),  # snr
            rng.uniform(0.0, 0.5, 60),  # swift_skew
            rng.uniform(-0.5, 0.5, 60),  # asym
            rng.uniform(20.0, 50.0, 60),  # phase_xy
            rng.uniform(-160.0, -120.0, 60),  # phase_yx
        ]
    )

    scorer = EMQCScorer(use_ml=True, n_estimators=10, random_state=0)
    assert scorer.fit(feat) is scorer
    assert scorer._model is not None

    scores = scorer.transform(feat)
    assert scores.shape == (60,)
    assert np.all((scores >= 0.0) & (scores <= 1.0))
    assert "ml+rules" in repr(scorer)


def test_emqc_scorer_score_table_and_site_transform(sites):
    from pycsamt.ai.processing.qc import EMQCScorer

    scorer = EMQCScorer(use_ml=False)
    scorer.fit(sites)

    tbl = scorer.score_table(sites)
    assert {"station", "freq", "score", "flag"}.issubset(tbl.columns)
    assert len(tbl) > 0
    assert set(np.unique(tbl["flag"])).issubset({0, 1})

    # transform() with a raw site collection (not ndarray/DataFrame) exercises
    # the `_to_feature_matrix` site-collection fallback branch.
    scores = scorer.transform(sites)
    assert scores.shape[0] == len(tbl)


def test_emqc_scorer_save_load_round_trip_with_ml_model(tmp_path):
    from pycsamt.ai.processing.qc import EMQCScorer

    rng = np.random.default_rng(2)
    feat = np.column_stack(
        [
            rng.uniform(2.0, 20.0, 40),
            rng.uniform(0.0, 0.5, 40),
            rng.uniform(-0.5, 0.5, 40),
            rng.uniform(20.0, 50.0, 40),
            rng.uniform(-160.0, -120.0, 40),
        ]
    )
    scorer = EMQCScorer(use_ml=True, n_estimators=10, random_state=0)
    scorer.fit(feat)
    scores = scorer.transform(feat)

    path = tmp_path / "qc_ml.npz"
    scorer.save(path)
    loaded = EMQCScorer.load(path)

    np.testing.assert_allclose(loaded.transform(feat), scores)


# ─────────────────────────────────────────────────────────────────────────────
# denoise.py
# ─────────────────────────────────────────────────────────────────────────────


def _denoise_X(n=24, n_comp=4, n_freqs=16, seed=5):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n, n_comp, n_freqs)).astype(np.float32)


def test_emdenoiser_verbose_training_prints(capsys):
    from pycsamt.ai.processing.denoise import EMDenoiser

    X = _denoise_X(n=16, n_freqs=12)
    den = EMDenoiser(channels=(4, 8, 4))
    den.fit(X, epochs=10, verbose=True)
    out = capsys.readouterr().out
    assert "Epoch" in out


def test_emdenoiser_save_load_round_trip_torch(tmp_path):
    from pycsamt.ai.processing.denoise import EMDenoiser

    X = _denoise_X(n=16, n_freqs=12)
    den = EMDenoiser(channels=(4, 8, 4))
    den.fit(X, epochs=2, verbose=False)
    assert den._backend_name == "torch"
    out_before = den.transform(X)

    path = tmp_path / "denoiser.npz"
    den.save(path)
    loaded = EMDenoiser.load(path)

    out_after = loaded.transform(X)
    np.testing.assert_allclose(out_before, out_after, rtol=1e-4, atol=1e-5)


def test_emdenoiser_numpy_fallback_scipy_missing(monkeypatch):
    import pycsamt.ai.processing.denoise as denoise_mod

    monkeypatch.setattr(denoise_mod, "active_backend", lambda: "none")
    monkeypatch.setitem(sys.modules, "scipy.ndimage", None)

    X = _denoise_X(n=8, n_freqs=10)
    den = denoise_mod.EMDenoiser()
    den.fit(X, epochs=1, verbose=False)
    assert den._use_numpy is True

    out = den.transform(X)
    assert out.shape == X.shape
    assert np.all(np.isfinite(out))


# ─────────────────────────────────────────────────────────────────────────────
# anomaly.py
# ─────────────────────────────────────────────────────────────────────────────


def _anomaly_X(n=40, n_feat=20, seed=6):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n, n_feat)).astype(np.float32)


def test_anomaly_detector_verbose_training_prints(capsys):
    from pycsamt.ai.processing.anomaly import AnomalyDetector

    X = _anomaly_X()
    det = AnomalyDetector(latent_dim=4, channels=(8,))
    det.fit(X, epochs=10, verbose=True)
    out = capsys.readouterr().out
    assert "AnomalyDetector" in out


def test_anomaly_detector_save_load_round_trip_torch(tmp_path):
    from pycsamt.ai.processing.anomaly import AnomalyDetector

    X = _anomaly_X()
    det = AnomalyDetector(latent_dim=4, channels=(8,))
    det.fit(X, epochs=2, verbose=False)
    assert det._backend_name == "torch"
    scores_before = det.transform(X)

    path = tmp_path / "anomaly_torch.npz"
    det.save(path)
    loaded = AnomalyDetector.load(path)

    scores_after = loaded.transform(X)
    np.testing.assert_allclose(scores_before, scores_after, rtol=1e-4, atol=1e-5)
    np.testing.assert_array_equal(loaded.flag_anomalies(X), det.flag_anomalies(X))


def test_anomaly_detector_pca_fallback_verbose_prints(monkeypatch, capsys):
    pytest.importorskip("sklearn")
    import pycsamt.ai.processing.anomaly as anomaly_mod

    monkeypatch.setattr(
        anomaly_mod.AnomalyDetector,
        "_fit_torch",
        lambda self, Xn, **kwargs: (_ for _ in ()).throw(ImportError("no torch")),
    )
    X = _anomaly_X(n=20, n_feat=10)
    det = anomaly_mod.AnomalyDetector(latent_dim=3)
    det.fit(X, epochs=1, verbose=True)
    out = capsys.readouterr().out
    assert "PCA fallback" in out
    assert det._use_pca is True


# ─────────────────────────────────────────────────────────────────────────────
# classify.py
# ─────────────────────────────────────────────────────────────────────────────


def _dim_Xy(n=80, seed=9):
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, 5)).astype(np.float32)
    X[:, 0] = np.abs(X[:, 0]) * 4
    X[:, 1] = np.abs(X[:, 1]) * 0.3
    from pycsamt.ai.processing.classify import _rule_labels

    y = _rule_labels(X[:, 0], X[:, 1])
    strike = np.where(y == 1, rng.uniform(-90.0, 90.0, n), np.nan).astype(np.float32)
    return X, y, strike


def test_dim_classifier_strike_training_torch_verbose(capsys):
    from pycsamt.ai.processing.classify import DimensionalityClassifier

    X, y, strike = _dim_Xy()
    clf = DimensionalityClassifier(hidden=(16,), dropout=0.0)
    clf.fit(X, y, strike=strike, epochs=10, verbose=True)
    out = capsys.readouterr().out
    assert "DimClassifier" in out
    assert clf._backend_name == "torch"

    pred_strike = clf.predict_strike(X)
    assert pred_strike.shape == (len(X),)
    labels = clf.predict(X)
    assert np.all(np.isnan(pred_strike[labels != 1]))


def test_dim_classifier_save_load_round_trip_torch(tmp_path):
    from pycsamt.ai.processing.classify import DimensionalityClassifier

    X, y, strike = _dim_Xy(n=60)
    clf = DimensionalityClassifier(hidden=(16,), dropout=0.0)
    clf.fit(X, y, strike=strike, epochs=2, verbose=False)
    proba_before = clf.transform(X)

    path = tmp_path / "dimclf_torch.npz"
    clf.save(path)
    loaded = DimensionalityClassifier.load(path)

    proba_after = loaded.transform(X)
    np.testing.assert_allclose(proba_before, proba_after, rtol=1e-4, atol=1e-5)
    np.testing.assert_array_equal(loaded.predict(X), clf.predict(X))


def test_dim_classifier_save_load_round_trip_rf(tmp_path, monkeypatch):
    pytest.importorskip("sklearn")
    import pycsamt.ai.processing.classify as classify_mod

    monkeypatch.setattr(
        classify_mod.DimensionalityClassifier,
        "_fit_torch",
        lambda self, Xn, y, strike, **kwargs: (_ for _ in ()).throw(
            ImportError("no torch")
        ),
    )
    X, y, _ = _dim_Xy(n=60)
    clf = classify_mod.DimensionalityClassifier()
    clf.fit(X, y, verbose=False)
    assert clf._use_rf is True
    proba_before = clf.transform(X)

    path = tmp_path / "dimclf_rf.npz"
    clf.save(path)
    loaded = classify_mod.DimensionalityClassifier.load(path)

    assert loaded._use_rf is True
    proba_after = loaded.transform(X)
    np.testing.assert_allclose(proba_before, proba_after)
    assert "rf" in repr(loaded)


def test_dim_classifier_predict_table_from_sites(sites):
    from pycsamt.ai.processing.classify import DimensionalityClassifier

    X, y, _ = _dim_Xy(n=60)
    clf = DimensionalityClassifier(hidden=(16,), dropout=0.0)
    clf.fit(X, y, epochs=2, verbose=False)

    tbl = clf.predict_table(sites)
    assert {
        "station",
        "freq",
        "period",
        "dim",
        "dim_label",
        "strike",
        "confidence",
    }.issubset(tbl.columns)
    assert len(tbl) > 0
    assert set(tbl["dim_label"].unique()).issubset({"1D", "2D", "3D"})


# ─────────────────────────────────────────────────────────────────────────────
# processing/plot.py — plot_qc_feature_heatmap, plot_qc_summary
# ─────────────────────────────────────────────────────────────────────────────


def _qc_full_dataframe():
    stations = [f"S{i:02d}" for i in range(6)]
    freqs = np.array([1000.0, 100.0, 10.0, 1.0])
    rng = np.random.default_rng(3)
    rows = []
    for st in stations:
        for fr in freqs:
            rows.append(
                {
                    "station": st,
                    "freq": fr,
                    "snr": rng.uniform(2, 30),
                    "swift_skew": rng.uniform(0, 0.5),
                    "asym": rng.uniform(-0.5, 0.5),
                    "phase_xy": rng.uniform(20, 50),
                    "phase_yx": rng.uniform(-160, -120),
                    "score": rng.uniform(0, 1),
                }
            )
    return pd.DataFrame(rows)


def test_plot_qc_feature_heatmap_default_and_overrides():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pycsamt.ai.processing.plot import plot_qc_feature_heatmap

    df = _qc_full_dataframe()

    fig = plot_qc_feature_heatmap(df, title="QC features")
    assert fig.__class__.__name__ == "Figure"
    assert len(fig.axes) >= 5
    plt.close(fig)

    fig2 = plot_qc_feature_heatmap(
        df,
        features=["snr", "asym"],
        cmaps={"snr": "viridis"},
        n_yticks=3,
    )
    assert len(fig2.axes) >= 2
    plt.close(fig2)

    # single-feature branch (n_feat == 1 -> axes wrapped in a list)
    fig3 = plot_qc_feature_heatmap(df, features=["snr"])
    plt.close(fig3)

    with pytest.raises(ValueError, match="No valid feature columns"):
        plot_qc_feature_heatmap(df, features=["not_a_column"])


def test_plot_qc_summary_dict_and_dataframe_inputs():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pycsamt.ai.processing.plot import plot_qc_summary

    scores = {
        "L18": np.array([0.2, 0.6, 0.9, 0.4]),
        "L22": np.array([0.7, 0.8, 0.3]),
    }
    fig = plot_qc_summary(scores, suptitle="QC summary", spread_kind="box")
    assert fig.__class__.__name__ == "Figure"
    assert len(fig.axes) >= 3
    plt.close(fig)

    df = _qc_full_dataframe()
    fig2 = plot_qc_summary(df, spread_kind="strip")
    plt.close(fig2)


if __name__ == "__main__":
    pytest.main([__file__])
