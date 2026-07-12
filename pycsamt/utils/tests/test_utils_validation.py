# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.utils.validation`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pycsamt.exceptions import NotFittedError, NotReadError
from pycsamt.utils.validation import (
    _assert_all_types,
    _is_arraylike_1d,
    _is_numeric_dtype,
    _validate_name_in,
    assert_ratio,
    check_consistency_size,
    check_has_read,
    check_is_fitted,
    ensure_n_items,
    has_read,
    isin,
    isin_if,
    isinstance_relaxed,
)


# ----------------------------- ensure_n_items -------------------------

def test_ensure_n_items_positional_and_iterable():
    assert ensure_n_items(1, 2) == (1, 2)
    assert ensure_n_items(items=[3, 4]) == (3, 4)
    assert ensure_n_items(items=(5, 6), return_as="list") == [5, 6]
    arr = ensure_n_items(items=[7, 8], return_as="array", dtype=float)
    assert isinstance(arr, np.ndarray)
    assert arr.dtype == float


def test_ensure_n_items_wrong_count_raises_or_warns():
    with pytest.raises(ValueError):
        ensure_n_items(1, 2, 3, n=2)
    with pytest.warns(UserWarning):
        out = ensure_n_items(1, 2, 3, n=2, error="warn")
    assert out == (1, 2)
    with pytest.warns(UserWarning):
        padded = ensure_n_items(items=[1], n=2, error="warn",
                                allow_none=True)
    assert padded == (1, None)


def test_ensure_n_items_numeric_policy():
    assert ensure_n_items("1", "2.5", expect="numeric",
                          coerce=True) == (1.0, 2.5)
    with pytest.raises(TypeError):
        ensure_n_items("a", "b", expect="numeric")
    with pytest.raises(ValueError):
        ensure_n_items(0.5, 9.0, expect="numeric", bounds=(0, 1))
    assert ensure_n_items(
        0.2, 0.8, expect="numeric", bounds=(0, 1)
    ) == (0.2, 0.8)


def test_ensure_n_items_string_policy_and_none():
    assert ensure_n_items("a", "b", expect="string") == ("a", "b")
    with pytest.raises(TypeError):
        ensure_n_items("a", 3, expect="string")
    with pytest.raises(TypeError):
        ensure_n_items(None, 1)
    assert ensure_n_items(None, 1, allow_none=True) == (None, 1)


def test_ensure_n_items_unique_and_to_string():
    with pytest.raises(ValueError):
        ensure_n_items(1, 1, unique=True)
    assert ensure_n_items(1, "x", to_string=True) == ("1", "x")


# ------------------------ has_read / check_has_read -------------------

class _Reader:
    def __init__(self):
        self.data = None

    def load(self):
        self.data = [1, 2, 3]

    @check_has_read(attributes="data")
    def total(self):
        return sum(self.data)


def test_has_read_attribute_inspection():
    r = _Reader()
    with pytest.raises(NotReadError):
        has_read(r, attributes="data")
    r.load()
    assert has_read(r, attributes="data") is True


def test_has_read_flag_and_custom_dunder():
    class Flagged:
        _has_read = False

    with pytest.raises(NotReadError):
        has_read(Flagged())
    Flagged._has_read = True
    assert has_read(Flagged()) is True

    class Dunder:
        def __has_read__(self):
            return False

    with pytest.raises(NotReadError):
        has_read(Dunder())


def test_has_read_requires_flag_or_attributes():
    class Empty:
        pass

    with pytest.raises(NotReadError):
        has_read(Empty())


def test_check_has_read_decorator_blocks_until_load():
    r = _Reader()
    with pytest.raises(NotReadError):
        r.total()
    r.load()
    assert r.total() == 6


# ----------------------------- check_is_fitted ------------------------

def test_check_is_fitted_variants():
    class WithDunder:
        def __is_fitted__(self):
            return True

    assert check_is_fitted(WithDunder()) is True

    class WithAttrs:
        coef_ = None

    with pytest.raises(NotFittedError):
        check_is_fitted(WithAttrs(), attributes=["coef_"])
    inst = WithAttrs()
    inst.coef_ = [1.0]
    assert check_is_fitted(inst, attributes=["coef_"]) is True

    class WithFlag:
        _is_fitted = False

    with pytest.raises(NotFittedError):
        check_is_fitted(WithFlag())

    class NoIndicator:
        pass

    with pytest.raises(NotFittedError):
        check_is_fitted(NoIndicator())


# ---------------------------- _assert_all_types -----------------------

def test_assert_all_types_pass_and_fail():
    assert _assert_all_types(3, int) == 3
    assert _assert_all_types(3.0, (int, float)) == 3.0
    with pytest.raises(TypeError):
        _assert_all_types("x", int, float, objname="param")
    with pytest.raises(TypeError):
        _assert_all_types(3)  # no expected types


# --------------------------------- isin -------------------------------

def test_isin_reduction_modes():
    assert isin([1, 2, 3], [2, 5], match="any") is True
    assert isin([1, 2, 3], [2, 5], match="all") is False
    assert isin([1, 2, 3], [2, 2, 3], match="count") == 2
    with pytest.raises(ValueError):
        isin([1], [1], match="most")


def test_isin_mask_and_extras():
    mask = isin([1, 2, 3], [2], return_mask=True)
    assert mask.tolist() == [False, True, False]
    inv = isin([1, 2, 3], [2], return_mask=True, invert=True)
    assert inv.tolist() == [True, False, True]


def test_isin_equal_nan_handling():
    # regression: equal_nan used to be forwarded to np.isin,
    # which does not accept it (TypeError on every call)
    assert isin([1.0, np.nan], [np.nan], match="all",
                equal_nan=True) is True
    assert isin([1.0, np.nan], [np.nan], match="all",
                equal_nan=False) is False
    mask = isin([1.0, np.nan], [np.nan], return_mask=True,
                equal_nan=True)
    assert mask.tolist() == [False, True]
    inv = isin([1.0, np.nan], [np.nan], return_mask=True,
               equal_nan=True, invert=True)
    assert inv.tolist() == [True, False]

    result, missing = isin([1, 2, 3], [2, 5], return_missing=True)
    assert result is False and missing == [5]
    result, missing, present = isin(
        [1, 2, 3], [2, 5], return_missing=True, return_present=True
    )
    assert present == [2]


# -------------------------------- isin_if ------------------------------

def test_isin_if_error_modes_and_returns():
    with pytest.raises(ValueError):
        isin_if(["a", "b"], ["c"])
    with pytest.warns(UserWarning):
        assert isin_if(["a", "b"], ["c"], error="warn") is None
    assert isin_if(["a", "b"], "a") is None
    assert isin_if(["a", "b"], ["a", "c"],
                   return_diff=True) == ["c"]
    assert isin_if(["a", "b"], ["a", "c"],
                   return_intersect=True) == ["a"]
    with pytest.raises(TypeError):
        isin_if(42, ["a"])


# ------------------------------ assert_ratio ---------------------------

def test_assert_ratio_values_and_percent():
    assert assert_ratio(0.5) == 0.5
    assert assert_ratio("25%") == pytest.approx(0.25)
    assert assert_ratio(30, in_percent=True) == pytest.approx(0.30)
    assert assert_ratio(0.3, bounds=(0, 1)) == pytest.approx(0.3)


def test_assert_ratio_errors():
    with pytest.raises(TypeError):
        assert_ratio("not-a-number")
    with pytest.raises(ValueError):
        assert_ratio(2.0, bounds=(0, 1))
    with pytest.raises(ValueError):
        assert_ratio(0.0, exclude_value=0.0)
    with pytest.raises(ValueError):
        assert_ratio(0.5, bounds=(0, 1, 2))


# ----------------------------- _validate_name_in ----------------------

def test_validate_name_in_modes():
    assert _validate_name_in("east", ("east", "north")) is True
    assert _validate_name_in("EAST ", ("east",),
                             expect_name="easting") == "easting"
    assert _validate_name_in("no", ("east",)) is False
    assert _validate_name_in("ast", "east", deep=True) is True
    with pytest.raises(KeyError):
        _validate_name_in("no", ("east",),
                          exception=KeyError("bad name"))


# ---------------------------- isinstance_relaxed ----------------------

def test_isinstance_relaxed_direct_and_by_name():
    class A:
        pass

    a = A()
    assert isinstance_relaxed(a, A)
    assert isinstance_relaxed(a, (int, A))
    assert not isinstance_relaxed(a, int)

    # simulate a module reload: same class name + module tail
    class B:
        pass

    B2 = type("B", (), {})
    B2.__module__ = B.__module__
    assert isinstance_relaxed(B2(), B)


# ------------------------- dtype / size helpers -----------------------

def test_is_numeric_dtype():
    assert _is_numeric_dtype(np.arange(3))
    assert _is_numeric_dtype([1, 2, 3], to_array=True)
    assert not _is_numeric_dtype(np.array(["a", "b"]))
    assert _is_numeric_dtype(pd.Series([1.0, 2.0]))
    with pytest.raises(TypeError):
        _is_numeric_dtype(42)
    with pytest.raises(ValueError):
        _is_numeric_dtype([1, 2])  # iterable but not array-like


def test_check_consistency_size():
    check_consistency_size(np.ones(3), [1, 2, 3], None)
    with pytest.raises(ValueError):
        check_consistency_size(np.ones(3), np.ones(4))


def test_is_arraylike_1d():
    assert _is_arraylike_1d(np.arange(3))
    assert _is_arraylike_1d(np.arange(3).reshape(3, 1))
    assert not _is_arraylike_1d(np.ones((2, 2)))
    # 0-d arrays intentionally count as 1-D (ndim < 2): callers
    # rely on this for scalar inputs wrapped by is_iterable
    assert _is_arraylike_1d(np.array(3.0))
    with pytest.raises(TypeError):
        _is_arraylike_1d("string")
