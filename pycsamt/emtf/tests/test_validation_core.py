from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtf.validation import (
    normalize_channels,
    normalize_periods,
    normalize_tf_data,
)


# ---------------------------------------------------------------------
# normalize_channels
# ---------------------------------------------------------------------


def test_normalize_channels_none_returns_empty_tuple():
    assert normalize_channels(None) == ()


def test_normalize_channels_strips_and_preserves_order():
    assert normalize_channels([" Hx ", "Hy"]) == ("Hx", "Hy")


def test_normalize_channels_rejects_empty_name():
    with pytest.raises(ValueError, match="non-empty"):
        normalize_channels(["Hx", "  "])


def test_normalize_channels_rejects_duplicates():
    with pytest.raises(ValueError, match="duplicate channel"):
        normalize_channels(["Hx", "Hx"])


# ---------------------------------------------------------------------
# normalize_periods
# ---------------------------------------------------------------------


def test_normalize_periods_none_returns_none():
    assert normalize_periods(None) is None


def test_normalize_periods_scalar_is_reshaped():
    out = normalize_periods(1.5)
    assert out.shape == (1,)
    assert out[0] == 1.5


def test_normalize_periods_rejects_non_1d():
    with pytest.raises(ValueError, match="1-D array"):
        normalize_periods(np.ones((2, 2)))


def test_normalize_periods_rejects_nonpositive_values():
    with pytest.raises(ValueError, match="finite positive"):
        normalize_periods([1.0, 0.0, 2.0])


def test_normalize_periods_rejects_non_finite_values():
    with pytest.raises(ValueError, match="finite positive"):
        normalize_periods([1.0, np.nan])


def test_normalize_periods_empty_array_skips_positivity_check():
    out = normalize_periods(np.array([]))
    assert out.shape == (0,)


def test_normalize_periods_enforces_expected_count():
    with pytest.raises(ValueError, match="period count"):
        normalize_periods([1.0, 2.0], n_periods=3)
    out = normalize_periods([1.0, 2.0], n_periods=2)
    assert out.size == 2


# ---------------------------------------------------------------------
# normalize_tf_data
# ---------------------------------------------------------------------


def test_normalize_tf_data_rejects_non_numeric():
    with pytest.raises(TypeError, match="must be numeric"):
        normalize_tf_data(["a", "b"], n_output=1, n_input=1)


def test_normalize_tf_data_scalar_promotes_for_1x1():
    out = normalize_tf_data(3.0, n_output=1, n_input=1)
    assert out.shape == (1, 1, 1)


def test_normalize_tf_data_scalar_rejected_for_non_1x1():
    with pytest.raises(ValueError, match="scalar 1x1"):
        normalize_tf_data(3.0, n_output=2, n_input=2)


def test_normalize_tf_data_1d_promotes_for_1x1():
    out = normalize_tf_data([1.0, 2.0, 3.0], n_output=1, n_input=1)
    assert out.shape == (3, 1, 1)


def test_normalize_tf_data_1d_rejected_for_non_1x1():
    with pytest.raises(ValueError, match="scalar input/output"):
        normalize_tf_data([1.0, 2.0], n_output=2, n_input=1)


def test_normalize_tf_data_2d_promotes_to_single_period():
    out = normalize_tf_data(np.ones((2, 2)), n_output=2, n_input=2)
    assert out.shape == (1, 2, 2)


def test_normalize_tf_data_2d_shape_mismatch_raises():
    with pytest.raises(ValueError, match="one TF matrix"):
        normalize_tf_data(np.ones((2, 3)), n_output=2, n_input=2)


def test_normalize_tf_data_3d_matches_channel_dims():
    out = normalize_tf_data(np.ones((5, 2, 2)), n_output=2, n_input=2)
    assert out.shape == (5, 2, 2)


def test_normalize_tf_data_3d_shape_mismatch_raises():
    with pytest.raises(ValueError, match="channel dimensions"):
        normalize_tf_data(np.ones((5, 1, 2)), n_output=2, n_input=2)


def test_normalize_tf_data_rejects_4d():
    with pytest.raises(ValueError, match="scalar, 1-D, 2-D, or 3-D"):
        normalize_tf_data(np.ones((2, 2, 2, 2)), n_output=2, n_input=2)
