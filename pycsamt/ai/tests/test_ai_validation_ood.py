"""Contracts for :mod:`pycsamt.ai.validation.ood`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.validation import (
    OODReport,
    flag_out_of_distribution,
    ood_score,
)

_REFERENCE = np.array(
    [
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [-1.0, 0.0],
        [0.0, -1.0],
        [0.5, 0.5],
    ]
)
_X = np.array([[0.0, 0.0], [50.0, 50.0]])


def test_mahalanobis_ranks_the_outlier_higher():
    scores = ood_score(_X, _REFERENCE, method="mahalanobis")
    assert scores.shape == (2,)
    assert scores[0] < scores[1]


def test_knn_ranks_the_outlier_higher():
    scores = ood_score(_X, _REFERENCE, method="knn", k=2)
    assert scores[0] < scores[1]


def test_ood_score_rejects_bad_method_and_mismatched_features():
    with pytest.raises(ValueError, match="method must be one of"):
        ood_score(_X, _REFERENCE, method="euclidean")
    with pytest.raises(ValueError, match="same number of features"):
        ood_score(np.zeros((2, 3)), _REFERENCE)


def test_validate_2d_rejects_wrong_ndim_empty_and_nonfinite():
    with pytest.raises(ValueError, match="2-D"):
        ood_score(np.zeros(3), _REFERENCE)
    with pytest.raises(ValueError, match="not be empty"):
        ood_score(np.zeros((0, 2)), _REFERENCE)
    bad = _X.copy()
    bad[0, 0] = np.nan
    with pytest.raises(ValueError, match="must be finite"):
        ood_score(bad, _REFERENCE)


def test_mahalanobis_requires_more_samples_than_features():
    small_reference = np.array([[0.0, 0.0], [1.0, 1.0]])
    with pytest.raises(ValueError, match="more samples than features"):
        ood_score(_X, small_reference, method="mahalanobis")


def test_mahalanobis_singular_covariance_raises():
    duplicate_reference = np.tile(np.array([[1.0, 2.0]]), (5, 1))
    with pytest.raises(ValueError, match="singular"):
        ood_score(_X, duplicate_reference, method="mahalanobis")


def test_knn_k_out_of_range_raises():
    with pytest.raises(ValueError, match="k must be an integer"):
        ood_score(_X, _REFERENCE, method="knn", k=0)
    with pytest.raises(ValueError, match="k must be an integer"):
        ood_score(_X, _REFERENCE, method="knn", k=7)


def test_flag_out_of_distribution_quantile_threshold():
    report = flag_out_of_distribution(
        _X, _REFERENCE, method="knn", k=2, quantile=0.5
    )
    assert isinstance(report, OODReport)
    assert report.quantile == pytest.approx(0.5)
    assert report.flagged.tolist() == [False, True]
    assert report.n_reference == 6
    assert report.n_features == 2
    assert report.fraction_flagged == pytest.approx(0.5)


def test_flag_out_of_distribution_explicit_threshold_overrides_quantile():
    report = flag_out_of_distribution(
        _X, _REFERENCE, method="knn", k=2, threshold=1000.0
    )
    assert report.quantile is None
    assert report.threshold == pytest.approx(1000.0)
    assert not report.flagged.any()


def test_flag_out_of_distribution_bad_quantile_and_threshold_raise():
    with pytest.raises(ValueError, match=r"within \(0, 1\)"):
        flag_out_of_distribution(_X, _REFERENCE, quantile=1.5)
    with pytest.raises(ValueError, match="threshold must be finite"):
        flag_out_of_distribution(_X, _REFERENCE, threshold=np.nan)


def test_flag_out_of_distribution_self_knn_k_out_of_range_raises():
    with pytest.raises(ValueError, match="self-referential"):
        flag_out_of_distribution(_X, _REFERENCE, method="knn", k=6)


def test_ood_report_validates_fields():
    with pytest.raises(ValueError, match="method must be one of"):
        OODReport(
            scores=np.array([1.0]),
            threshold=1.0,
            flagged=np.array([False]),
            method="bogus",
            quantile=None,
            n_reference=5,
            n_features=2,
        )
    with pytest.raises(ValueError, match="flagged must have"):
        OODReport(
            scores=np.array([1.0, 2.0]),
            threshold=1.0,
            flagged=np.array([False]),
            method="knn",
            quantile=None,
            n_reference=5,
            n_features=2,
        )
    with pytest.raises(ValueError, match="threshold must be finite"):
        OODReport(
            scores=np.array([1.0]),
            threshold=np.nan,
            flagged=np.array([False]),
            method="knn",
            quantile=None,
            n_reference=5,
            n_features=2,
        )
    with pytest.raises(ValueError, match=r"within \(0, 1\)"):
        OODReport(
            scores=np.array([1.0]),
            threshold=1.0,
            flagged=np.array([False]),
            method="knn",
            quantile=1.5,
            n_reference=5,
            n_features=2,
        )
    with pytest.raises(ValueError, match="must be positive"):
        OODReport(
            scores=np.array([1.0]),
            threshold=1.0,
            flagged=np.array([False]),
            method="knn",
            quantile=None,
            n_reference=0,
            n_features=2,
        )
