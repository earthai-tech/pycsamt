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
        padded = ensure_n_items(items=[1], n=2, error="warn", allow_none=True)
    assert padded == (1, None)


def test_ensure_n_items_numeric_policy():
    assert ensure_n_items("1", "2.5", expect="numeric", coerce=True) == (
        1.0,
        2.5,
    )
    with pytest.raises(TypeError):
        ensure_n_items("a", "b", expect="numeric")
    with pytest.raises(ValueError):
        ensure_n_items(0.5, 9.0, expect="numeric", bounds=(0, 1))
    assert ensure_n_items(0.2, 0.8, expect="numeric", bounds=(0, 1)) == (
        0.2,
        0.8,
    )


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


def test_ensure_n_items_scalar_items_is_wrapped_in_a_list():
    # items=None is not iterable -> wrapped as [None]
    out = ensure_n_items(items=None, n=1, allow_none=True)
    assert out == (None,)


def test_ensure_n_items_too_many_warn_truncates():
    with pytest.warns(UserWarning, match="truncating"):
        out = ensure_n_items(1, 2, 3, n=2, error="warn")
    assert out == (1, 2)


def test_ensure_n_items_none_warn_mode_keeps_none():
    with pytest.warns(UserWarning, match="not allowed"):
        out = ensure_n_items(None, 1, error="warn")
    assert out == (None, 1)


def test_ensure_n_items_string_policy_warn_mode_coerces_anyway():
    with pytest.warns(UserWarning, match="must be string"):
        out = ensure_n_items("a", 3, expect="string", error="warn")
    assert out == ("a", "3")


def test_ensure_n_items_numeric_policy_warn_mode_paths():
    with pytest.warns(UserWarning, match="not numeric"):
        out = ensure_n_items("a", "b", expect="numeric", error="warn")
    assert all(np.isnan(v) for v in out)

    with pytest.warns(UserWarning, match="is NaN"):
        out = ensure_n_items(
            float("nan"), 1.0, expect="numeric", error="warn",
        )
    assert np.isnan(out[0])

    with pytest.warns(UserWarning, match="outside"):
        out = ensure_n_items(
            0.5, 9.0, expect="numeric", bounds=(0, 1), error="warn",
        )
    assert out == (0.5, 9.0)


def test_ensure_n_items_unique_warn_mode_keeps_duplicates():
    with pytest.warns(UserWarning, match="must be unique"):
        out = ensure_n_items(1, 1, unique=True, error="warn")
    assert out == (1, 1)


def test_ensure_n_items_error_ignore_suppresses_all_diagnostics(recwarn):
    out = ensure_n_items(1, 2, 3, n=2, error="ignore")
    assert out == (1, 2)
    assert len(recwarn) == 0


def test_ensure_n_items_error_ignore_suppresses_every_element_check(recwarn):
    assert ensure_n_items(None, 1, error="ignore") == (None, 1)
    assert ensure_n_items(
        "a", 3, expect="string", error="ignore",
    ) == ("a", "3")
    out = ensure_n_items("a", "b", expect="numeric", error="ignore")
    assert all(np.isnan(v) for v in out)
    out = ensure_n_items(
        float("nan"), 1.0, expect="numeric", error="ignore",
    )
    assert np.isnan(out[0])
    out = ensure_n_items(
        0.5, 9.0, expect="numeric", bounds=(0, 1), error="ignore",
    )
    assert out == (0.5, 9.0)
    assert ensure_n_items(1, 1, unique=True, error="ignore") == (1, 1)
    assert len(recwarn) == 0


def test_ensure_n_items_unique_true_with_actually_unique_values():
    assert ensure_n_items(1, 2, unique=True) == (1, 2)


def test_ensure_n_items_string_items_is_treated_as_scalar():
    # A bare string is not split into characters.
    out = ensure_n_items(items="ab", n=1)
    assert out == ("ab",)


# ------------------------------- _to_float -----------------------------


def test_to_float_direct_branches():
    from pycsamt.utils.validation import _to_float

    assert np.isnan(_to_float(None))
    assert _to_float(3.5) == 3.5
    assert _to_float(np.float64(2.0)) == 2.0
    assert _to_float(4) == 4.0
    assert _to_float(np.int32(5)) == 5.0
    with pytest.raises(TypeError, match="bool is not accepted"):
        _to_float(True)
    with pytest.raises(TypeError, match="Cannot coerce"):
        _to_float("not-a-number")


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


def test_has_read_non_bool_flag_falls_through_to_attributes():
    class Reader:
        _has_read = "maybe"  # not a bool -> falls through
        data = [1, 2]

    assert has_read(Reader(), attributes="data") is True


def test_has_read_dunder_true_returns_true_directly():
    class Dunder:
        def __has_read__(self):
            return True

    assert has_read(Dunder()) is True


def test_has_read_finds_self_from_caller_frame_when_obj_omitted():
    class Reader:
        def __init__(self):
            self.data = [1]

        def check(self):
            return has_read(attributes="data")

    assert Reader().check() is True


def test_has_read_raises_when_no_object_and_no_caller_self():
    with pytest.raises(ValueError, match="No object provided"):
        has_read()


def test_has_read_empty_dataframe_attribute_is_not_read():
    class Reader:
        df = pd.DataFrame()

    with pytest.raises(NotReadError):
        has_read(Reader(), attributes="df")


def test_has_read_empty_sequence_attribute_is_not_read():
    class Reader:
        items = []

    with pytest.raises(NotReadError):
        has_read(Reader(), attributes="items")


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


def test_check_is_fitted_finds_self_from_caller_frame():
    class Model:
        _is_fitted = True

        def check(self):
            return check_is_fitted()

    assert Model().check() is True


def test_check_is_fitted_fitted_flag_fallback_name():
    class Model:
        _fitted = True

    assert check_is_fitted(Model()) is True

    class NotFitted:
        _fitted = False

    with pytest.raises(NotFittedError):
        check_is_fitted(NotFitted())

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
    assert isin([1.0, np.nan], [np.nan], match="all", equal_nan=True) is True
    assert isin([1.0, np.nan], [np.nan], match="all", equal_nan=False) is False
    mask = isin([1.0, np.nan], [np.nan], return_mask=True, equal_nan=True)
    assert mask.tolist() == [False, True]
    inv = isin([1.0, np.nan], [np.nan], return_mask=True, equal_nan=True, invert=True)
    assert inv.tolist() == [True, False]

    result, missing = isin([1, 2, 3], [2, 5], return_missing=True)
    assert result is False and missing == [5]
    result, missing, present = isin(
        [1, 2, 3], [2, 5], return_missing=True, return_present=True
    )
    assert present == [2]


def test_isin_count_mode_falls_back_when_unique_raises():
    # object-dtype arrays with mixed, mutually-incomparable types make
    # np.unique raise TypeError; the count path must still work via
    # dict.fromkeys-based de-duplication.
    result = isin([1, "a", 2], [1, "a", "a"], match="count")
    assert result == 2


# -------------------------------- _isin ---------------------------------


def test_private_isin_all_present_and_missing():
    from pycsamt.utils.validation import _isin

    assert _isin([1, 2, 3], [1, 2]) is True
    assert _isin([1, 2, 3], [5]) is False


def test_private_isin_return_mask():
    from pycsamt.utils.validation import _isin

    mask = _isin([1, 2, 3], [2], return_mask=True)
    assert mask.tolist() == [False, True, False]


def test_private_isin_scalar_subarr():
    from pycsamt.utils.validation import _isin

    assert _isin([1, 2, 3], 2) is True
    assert _isin([1, 2, 3], 9) is False


class _UnconvertibleToArray:
    """Raises when NumPy tries to coerce it, simulating a bad input."""

    def __array__(self, dtype=None):
        raise RuntimeError("cannot convert")


def test_private_isin_unique_failure_still_returns_correct_result():
    from pycsamt.utils.validation import _isin

    # mixed, mutually-incomparable types make np.unique raise internally;
    # _isin must fall back to the un-deduplicated array and still work.
    assert _isin([1, "a", 2], [1, "a", "a"]) is True


def test_private_isin_rejects_unconvertible_input():
    from pycsamt.utils.validation import _isin

    with pytest.raises(ValueError, match="Invalid inputs for membership"):
        _isin(_UnconvertibleToArray(), [1])


def test_isin_rejects_unconvertible_input():
    with pytest.raises(ValueError, match="Invalid inputs for membership"):
        isin(_UnconvertibleToArray(), [1])


# -------------------------------- isin_if ------------------------------


def test_isin_if_error_modes_and_returns():
    with pytest.raises(ValueError):
        isin_if(["a", "b"], ["c"])
    with pytest.warns(UserWarning):
        assert isin_if(["a", "b"], ["c"], error="warn") is None
    assert isin_if(["a", "b"], "a") is None
    assert isin_if(["a", "b"], ["a", "c"], return_diff=True) == ["c"]
    assert isin_if(["a", "b"], ["a", "c"], return_intersect=True) == ["a"]
    with pytest.raises(TypeError):
        isin_if(42, ["a"])


def test_isin_if_non_iterable_items_becomes_a_singleton_set():
    # items=5 is not iterable -> set(5) raises TypeError internally,
    # caught and treated as the scalar singleton {5}.
    assert isin_if([1, 2, 5], 5) is None  # no error: 5 is present
    with pytest.raises(ValueError):
        isin_if([1, 2, 3], 5)


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


def test_assert_ratio_percent_value_already_a_fraction_is_not_rescaled():
    # 0 < val <= 1 with in_percent=True does not get divided by 100.
    assert assert_ratio(0.5, in_percent=True) == pytest.approx(0.5)


def test_assert_ratio_bounds_with_one_open_end_skips_range_check():
    assert assert_ratio(500.0, bounds=(0, None)) == 500.0
    assert assert_ratio(-500.0, bounds=(None, 0)) == -500.0


def test_assert_ratio_exclude_value_conversion_failure_falls_back_to_low():
    with pytest.raises(ValueError, match="excluding"):
        assert_ratio(0.0, bounds=(0, 1), exclude_value="not-a-number")


def test_assert_ratio_exclude_value_failure_without_bounds_warns():
    with pytest.warns(UserWarning, match="Cannot exclude"):
        assert_ratio(0.5, exclude_value="not-a-number")


def test_assert_ratio_percent_result_over_one_raises():
    with pytest.raises(ValueError, match="must be <= 1.0"):
        assert_ratio(150, in_percent=True, bounds=(0, 200))


# ----------------------------- _validate_name_in ----------------------


def test_validate_name_in_modes():
    assert _validate_name_in("east", ("east", "north")) is True
    assert _validate_name_in("EAST ", ("east",), expect_name="easting") == "easting"
    assert _validate_name_in("no", ("east",)) is False
    assert _validate_name_in("ast", "east", deep=True) is True
    with pytest.raises(KeyError):
        _validate_name_in("no", ("east",), exception=KeyError("bad name"))


def test_validate_name_in_rejects_unconvertible_defaults():
    with pytest.raises(TypeError, match="must be str or sequence"):
        _validate_name_in("east", defaults=42)


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


def test_isinstance_relaxed_tuple_target_with_relaxed_match():
    class C:
        pass

    C2 = type("C", (), {})
    C2.__module__ = C.__module__
    # fast isinstance() fails (different class objects); relaxed
    # name+module-tail matching over a tuple target must still pass.
    assert isinstance_relaxed(C2(), (int, C))


def test_isinstance_relaxed_skips_unnamed_targets():
    class _NoNameMeta(type):
        @property
        def __name__(cls):  # noqa: N807 - intentional, simulates a
            # type object whose __name__ access raises (some C-extension
            # or proxy types behave this way).
            raise AttributeError("no name")

    class Unnamed(metaclass=_NoNameMeta):
        pass

    class D:
        pass

    d = D()
    # Unnamed's __name__ is unresolvable, so it must be skipped rather
    # than erroring.
    assert not isinstance_relaxed(d, (Unnamed,))


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


def test_private_check_consistency_size_matches_and_mismatches():
    from pycsamt.utils.validation import _check_consistency_size

    assert _check_consistency_size([1, 2, 3], [4, 5, 6]) is True
    with pytest.raises(AssertionError, match="Array sizes must match"):
        _check_consistency_size([1, 2], [1, 2, 3])
    assert _check_consistency_size([1, 2], [1, 2, 3], error="ignore") is False


def test_is_arraylike_1d():
    assert _is_arraylike_1d(np.arange(3))
    assert _is_arraylike_1d(np.arange(3).reshape(3, 1))
    assert not _is_arraylike_1d(np.ones((2, 2)))
    # 0-d arrays intentionally count as 1-D (ndim < 2): callers
    # rely on this for scalar inputs wrapped by is_iterable
    assert _is_arraylike_1d(np.array(3.0))
    with pytest.raises(TypeError):
        _is_arraylike_1d("string")


def test_is_arraylike_1d_numpy_scalar_is_not_1d():
    # np.float64 has __array__ but np.isscalar() is True for it,
    # so it must be rejected as a scalar rather than treated as 1-D.
    assert not _is_arraylike_1d(np.float64(3.0))


def test_is_arraylike_1d_falls_back_to_len_when_ndim_shape_missing():
    class ArrayLikeNoNdim:
        def __array__(self, dtype=None):
            return np.array([1, 2, 3])

        def __len__(self):
            return 3

    assert _is_arraylike_1d(ArrayLikeNoNdim())

    class ArrayLikeUnsizeable:
        def __array__(self, dtype=None):
            return np.array([1, 2, 3])

        def __len__(self):
            raise RuntimeError("no length")

    assert not _is_arraylike_1d(ArrayLikeUnsizeable())
