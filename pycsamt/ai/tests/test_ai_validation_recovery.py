"""Contracts for :mod:`pycsamt.ai.validation.recovery`."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.ai.validation import (
    RecoveryReport,
    depth_profile_mae,
    depth_profile_rmse,
    recovery_report,
    structural_similarity,
)

_PRED = np.array([[1.0, 2.0], [3.0, 4.0]])
_TRUE = np.array([[1.0, 2.0], [3.0, 6.0]])


def test_recovery_report_basic_values():
    report = recovery_report(_PRED, _TRUE)
    assert isinstance(report, RecoveryReport)
    assert report.rmse == pytest.approx(1.0)
    assert report.mae == pytest.approx(0.5)
    assert report.r2 == pytest.approx(1.0 - 4.0 / 14.0)
    assert report.n_valid == 4
    assert report.shape == (2, 2)
    np.testing.assert_allclose(
        report.depth_rmse, [0.0, np.sqrt(2.0)]
    )
    np.testing.assert_allclose(report.depth_mae, [0.0, 1.0])
    # grid is smaller than the default ssim_window=7: skipped.
    assert report.ssim is None


def test_recovery_report_computes_ssim_on_a_large_clean_grid():
    grid = np.arange(64, dtype=float).reshape(8, 8)
    report = recovery_report(grid, grid, ssim_window=3)
    assert report.ssim == pytest.approx(1.0)
    assert report.rmse == pytest.approx(0.0)


def test_recovery_report_skips_ssim_when_partially_masked():
    grid = np.arange(64, dtype=float).reshape(8, 8)
    valid = np.ones((8, 8), dtype=bool)
    valid[0, 0] = False
    report = recovery_report(grid, grid, valid=valid, ssim_window=3)
    assert report.ssim is None
    assert report.n_valid == 63


def test_recovery_report_r2_is_nan_for_constant_truth():
    pred = np.array([[1.0, 1.0], [1.0, 1.0]])
    true = np.array([[2.0, 2.0], [2.0, 2.0]])
    report = recovery_report(pred, true, compute_ssim=False)
    assert np.isnan(report.r2)


def test_recovery_report_shape_ndim_and_empty_checks_raise():
    with pytest.raises(ValueError, match="share one shape"):
        recovery_report(np.zeros((2, 2)), np.zeros((2, 3)))
    with pytest.raises(ValueError, match="2-D or 3-D"):
        recovery_report(np.zeros(4), np.zeros(4))
    with pytest.raises(ValueError, match="not be empty"):
        recovery_report(np.zeros((0, 2)), np.zeros((0, 2)))


def test_recovery_report_valid_mask_and_no_valid_cells():
    result = recovery_report(
        _PRED,
        _TRUE,
        valid=np.array([[True, True], [True, False]]),
        compute_ssim=False,
    )
    assert result.n_valid == 3
    assert result.rmse == pytest.approx(0.0)
    with pytest.raises(ValueError, match="no valid cell"):
        recovery_report(
            _PRED, _TRUE, valid=np.zeros((2, 2), dtype=bool)
        )


def test_depth_profile_rmse_and_mae_support_negative_and_bad_axis():
    rmse_axis0 = depth_profile_rmse(_PRED, _TRUE, axis=0)
    rmse_axis_neg = depth_profile_rmse(_PRED, _TRUE, axis=-2)
    np.testing.assert_allclose(rmse_axis0, rmse_axis_neg)
    mae = depth_profile_mae(_PRED, _TRUE, axis=1)
    assert mae.shape == (2,)
    with pytest.raises(ValueError, match="axis must be in"):
        depth_profile_rmse(_PRED, _TRUE, axis=5)


def test_depth_profile_layer_with_no_valid_cells_is_nan():
    valid = np.array([[True, True], [False, False]])
    profile = depth_profile_rmse(_PRED, _TRUE, axis=0, valid=valid)
    assert profile[0] == pytest.approx(0.0)
    assert np.isnan(profile[1])


def test_structural_similarity_identical_vs_shifted_grids():
    grid = np.arange(64, dtype=float).reshape(8, 8)
    same = structural_similarity(grid, grid, window=3)
    shifted = structural_similarity(grid, grid + 5.0, window=3)
    assert same == pytest.approx(1.0)
    assert shifted < same


def test_structural_similarity_rejects_bad_window_and_nonfinite():
    grid = np.arange(64, dtype=float).reshape(8, 8)
    with pytest.raises(ValueError, match="odd integer"):
        structural_similarity(grid, grid, window=4)
    with pytest.raises(ValueError, match="odd integer"):
        structural_similarity(grid, grid, window=0)
    with pytest.raises(ValueError, match="smallest grid axis"):
        structural_similarity(grid, grid, window=9)
    bad = grid.copy()
    bad[0, 0] = np.nan
    with pytest.raises(ValueError, match="fully finite"):
        structural_similarity(bad, grid, window=3)


def test_structural_similarity_rejects_bad_data_range():
    grid = np.arange(64, dtype=float).reshape(8, 8)
    with pytest.raises(ValueError, match="data_range"):
        structural_similarity(grid, grid, window=3, data_range=0.0)
