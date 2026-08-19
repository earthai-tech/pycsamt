from __future__ import annotations

import numpy as np
import pytest

from pycsamt.emtf import StatisticalEstimate


def test_construction_normalizes_name_and_kind():
    est = StatisticalEstimate(
        name=" var ",
        data=[1.0, 2.0, 3.0],
        kind=" Variance ",
    )
    assert est.name == "VAR"
    assert est.kind == "variance"
    assert est.units is None
    assert est.attrs == {}


def test_empty_name_raises():
    with pytest.raises(ValueError, match="name must be non-empty"):
        StatisticalEstimate(name="  ", data=[1.0], kind="variance")


def test_empty_kind_raises():
    with pytest.raises(ValueError, match="kind must be non-empty"):
        StatisticalEstimate(name="VAR", data=[1.0], kind="  ")


def test_non_numeric_data_raises_type_error():
    with pytest.raises(TypeError, match="must be numeric"):
        StatisticalEstimate(name="VAR", data=["a", "b"], kind="variance")


def test_attrs_none_becomes_empty_dict():
    est = StatisticalEstimate(
        name="VAR", data=[1.0], kind="variance", attrs=None
    )
    assert est.attrs == {}


def test_shape_and_is_complex_real():
    est = StatisticalEstimate(
        name="VAR", data=np.zeros((3, 2, 2)), kind="variance"
    )
    assert est.shape == (3, 2, 2)
    assert est.is_complex is False


def test_is_complex_true_for_complex_dtype():
    est = StatisticalEstimate(
        name="INVSIGCOV",
        data=np.zeros((2, 2), dtype=complex),
        kind="inverse_signal_covariance",
    )
    assert est.is_complex is True


def test_copy_is_detached_from_original():
    original = StatisticalEstimate(
        name="VAR", data=np.array([1.0, 2.0]), kind="variance"
    )
    duplicate = original.copy()
    duplicate.data[0] = 99.0
    assert original.data[0] == 1.0
    assert duplicate is not original
    assert duplicate.attrs is not original.attrs
