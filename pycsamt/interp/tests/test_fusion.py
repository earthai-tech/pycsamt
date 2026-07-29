# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.interp.fusion — MultiMethodEMModel."""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.interp._base import ResistivityModel
from pycsamt.interp.fusion import (
    FusionDiagnostics,
    MultiMethodEMModel,
    _interp_model_to_grid,
    _sigmoid,
)


def _primary(rms=2.0, station_x=None, station_names=None):
    """Shallow model: z in [10, 70], rho = 100 Ohm.m everywhere."""
    rho = np.full((4, 2), np.log10(100.0))
    return ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 100.0]),
        z_centers=np.array([10.0, 30.0, 50.0, 70.0]),
        station_x=station_x,
        station_names=station_names,
        method="TDEM",
        rms=rms,
    )


def _secondary(rms=4.0):
    """Deep model: z in [60, 300], rho = 1000 Ohm.m everywhere."""
    rho = np.full((4, 2), np.log10(1000.0))
    return ResistivityModel.from_array(
        rho,
        x_centers=np.array([0.0, 100.0]),
        z_centers=np.array([60.0, 120.0, 200.0, 300.0]),
        method="AMT",
        rms=rms,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Construction
# ─────────────────────────────────────────────────────────────────────────────


def test_invalid_blend_raises():
    with pytest.raises(ValueError, match="blend must be one of"):
        MultiMethodEMModel(_primary(), _secondary(), blend="bogus")


def test_repr():
    m = MultiMethodEMModel(_primary(), _secondary(), blend="sigmoid")
    r = repr(m)
    assert "TDEM" in r and "AMT" in r and "sigmoid" in r


def test_diagnostics_none_before_merge():
    m = MultiMethodEMModel(_primary(), _secondary())
    assert m.diagnostics_ is None


# ─────────────────────────────────────────────────────────────────────────────
# merge() — linear blend (default), with overlap
# ─────────────────────────────────────────────────────────────────────────────


def test_merge_linear_basic():
    m = MultiMethodEMModel(_primary(), _secondary(), blend="linear")
    fused = m.merge()
    assert fused.method == "TDEM+AMT"
    assert np.isnan(fused.rms)
    assert m.diagnostics_ is not None
    assert isinstance(m.diagnostics_, FusionDiagnostics)
    assert m.diagnostics_.has_overlap is True
    assert m.diagnostics_.blend_mode == "linear"
    assert m.diagnostics_.primary_method == "TDEM"
    assert m.diagnostics_.secondary_method == "AMT"
    assert m.diagnostics_.n_z_fused == len(fused.z_centers)
    # weights: 1 at shallow, 0 at deep
    assert m.diagnostics_.blend_weights[0] == pytest.approx(1.0)
    assert m.diagnostics_.blend_weights[-1] == pytest.approx(0.0)


def test_merge_shallowest_cell_is_pure_primary():
    m = MultiMethodEMModel(_primary(), _secondary(), blend="linear")
    fused = m.merge()
    # z=10 is well above the overlap zone [60,70] -> pure primary (100 Ohm.m)
    iz = list(fused.z_centers).index(10.0)
    assert fused.rho_2d[iz, 0] == pytest.approx(np.log10(100.0))


def test_merge_deepest_cell_is_pure_secondary():
    m = MultiMethodEMModel(_primary(), _secondary(), blend="linear")
    fused = m.merge()
    iz = list(fused.z_centers).index(300.0)
    assert fused.rho_2d[iz, 0] == pytest.approx(np.log10(1000.0))


def test_merge_station_metadata_from_primary_when_present():
    p = _primary(station_x=np.array([0.0, 100.0]), station_names=["A", "B"])
    m = MultiMethodEMModel(p, _secondary())
    fused = m.merge()
    assert fused.station_names == ["A", "B"]
    np.testing.assert_array_equal(fused.station_x, [0.0, 100.0])


def test_merge_station_metadata_falls_back_to_x_out_when_empty():
    p = _primary(station_x=np.array([]), station_names=[])
    m = MultiMethodEMModel(p, _secondary())
    fused = m.merge()
    np.testing.assert_array_equal(fused.station_x, fused.x_centers)


# ─────────────────────────────────────────────────────────────────────────────
# merge() — sigmoid blend
# ─────────────────────────────────────────────────────────────────────────────


def test_merge_sigmoid_blend_monotonic_weights():
    m = MultiMethodEMModel(_primary(), _secondary(), blend="sigmoid")
    m.merge()
    w = m.diagnostics_.blend_weights
    # weights must be non-increasing with depth
    assert np.all(np.diff(w) <= 1e-9)


# ─────────────────────────────────────────────────────────────────────────────
# merge() — rms_weighted blend
# ─────────────────────────────────────────────────────────────────────────────


def test_merge_rms_weighted_constant_in_overlap():
    p = _primary(rms=2.0)
    s = _secondary(rms=4.0)
    m = MultiMethodEMModel(p, s, blend="rms_weighted")
    m.merge()
    w = m.diagnostics_.blend_weights
    expected_w = 4.0 / (2.0 + 4.0)
    # cells strictly inside the overlap [60,70) should carry the constant weight
    z_out = np.array(
        sorted(
            set(p.z_centers[p.z_centers <= 70.0])
            | set(s.z_centers[(s.z_centers >= 60.0) & (s.z_centers <= 70.0)])
            | set(s.z_centers[s.z_centers > 70.0])
        )
    )
    in_overlap = (z_out > 60.0) & (z_out < 70.0)
    if in_overlap.any():
        assert w[in_overlap][0] == pytest.approx(expected_w)


def test_merge_rms_weighted_falls_back_to_linear_when_nan():
    p = _primary(rms=float("nan"))
    s = _secondary(rms=4.0)
    m = MultiMethodEMModel(p, s, blend="rms_weighted")
    m.merge()
    # should not crash, and behave like linear (weights vary across overlap)
    w = m.diagnostics_.blend_weights
    assert w[0] == pytest.approx(1.0)
    assert w[-1] == pytest.approx(0.0)


def test_merge_rms_weighted_falls_back_when_sum_zero():
    p = _primary(rms=0.0)
    s = _secondary(rms=0.0)
    m = MultiMethodEMModel(p, s, blend="rms_weighted")
    m.merge()  # must not raise ZeroDivisionError
    w = m.diagnostics_.blend_weights
    assert w[0] == pytest.approx(1.0)


# ─────────────────────────────────────────────────────────────────────────────
# merge() — no overlap
# ─────────────────────────────────────────────────────────────────────────────


def test_merge_no_overlap_hard_boundary():
    p = _primary()
    s = _secondary()
    m = MultiMethodEMModel(p, s, primary_max_depth=50.0, secondary_min_depth=100.0)
    fused = m.merge()
    assert m.diagnostics_.has_overlap is False
    w = m.diagnostics_.blend_weights
    z = fused.z_centers
    # hard boundary is at z_lo = secondary_min_depth (z_ov_start): primary
    # weight is 1 at/below it, 0 above it
    assert np.all(w[z <= 100.0] == 1.0)
    assert np.all(w[z > 100.0] == 0.0)


# ─────────────────────────────────────────────────────────────────────────────
# blend_overlap / z_grid / depth overrides
# ─────────────────────────────────────────────────────────────────────────────


def test_merge_blend_overlap_narrows_window():
    m = MultiMethodEMModel(_primary(), _secondary(), blend_overlap=2.0)
    m.merge()
    assert (
        m.diagnostics_.z_overlap_end - m.diagnostics_.z_overlap_start
        == pytest.approx(2.0)
    )


def test_merge_explicit_z_grid():
    z_grid = np.array([0.0, 50.0, 100.0, 250.0])
    m = MultiMethodEMModel(_primary(), _secondary(), z_grid=z_grid)
    fused = m.merge()
    np.testing.assert_array_equal(fused.z_centers, z_grid)


def test_merge_primary_max_depth_override_changes_grid():
    m1 = MultiMethodEMModel(_primary(), _secondary())
    m2 = MultiMethodEMModel(_primary(), _secondary(), primary_max_depth=30.0)
    fused1 = m1.merge()
    fused2 = m2.merge()
    assert not np.array_equal(fused1.z_centers, fused2.z_centers)


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────


def test_interp_model_to_grid_shape_and_values():
    p = _primary()
    z_out = np.array([10.0, 40.0, 70.0])
    x_out = np.array([0.0, 50.0, 100.0])
    rho = _interp_model_to_grid(p, x_out, z_out)
    assert rho.shape == (3, 3)
    np.testing.assert_allclose(rho, np.log10(100.0))


def test_sigmoid_midpoint_and_bounds():
    t = np.array([0.0, 0.5, 1.0])
    w = _sigmoid(t, k=10.0)
    assert w[1] == pytest.approx(0.5)
    assert w[0] < 0.5 < w[2]


def test_sigmoid_monotonic():
    t = np.linspace(0, 1, 20)
    w = _sigmoid(t, k=5.0)
    assert np.all(np.diff(w) >= 0)
