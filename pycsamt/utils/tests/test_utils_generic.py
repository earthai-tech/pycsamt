# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.utils.generic`.

``ensure_package`` is not exercised here: it shells out to pip and
must never run inside the offline test suite.
"""

from __future__ import annotations

import numpy as np
import pytest

from pycsamt.utils.generic import (
    count_functions,
    error_policy,
    get_valid_kwargs,
    make_ids,
    strip_item,
)


# ------------------------------- make_ids ------------------------------

def test_make_ids_from_sequence_and_count():
    assert make_ids(range(3), prefix="ix") == ["ix0", "ix1", "ix2"]
    # width follows the number of digits of the largest index
    assert make_ids(8, prefix="L", start=1)[:3] == ["L1", "L2", "L3"]
    assert make_ids(10, prefix="L", start=1)[0] == "L01"
    assert make_ids(10, prefix="L", start=1)[-1] == "L10"


def test_make_ids_zfill_and_skip_leading_zero():
    ids = make_ids(12, prefix="line", zfill=4, start=1)
    assert ids[0] == "line0001"
    assert ids[-1] == "line0012"
    assert make_ids(5, prefix="S", skip_leading_zero=True) == [
        "S0", "S1", "S2", "S3", "S4"
    ]


def test_make_ids_numpy_count_and_negative():
    assert make_ids(np.int64(2), prefix="n") == ["n0", "n1"]
    with pytest.raises(ValueError):
        make_ids(-1)


# ------------------------------ strip_item -----------------------------

def test_strip_item_whitespace_paths():
    assert strip_item("    ss_data   ") == "ss_data"
    assert strip_item(["  a  ", "   b"]) == ["a", "b"]
    assert strip_item(None) is None


def test_strip_item_token_and_array():
    arr = np.array(["////name////", "////x"], dtype="<U16")
    out = strip_item(arr, item="//")
    assert isinstance(out, np.ndarray)
    assert out.tolist() == ["name", "x"]


def test_strip_item_all_blank_returns_none():
    with pytest.warns(RuntimeWarning):
        assert strip_item(["      ", "   "]) is None


def test_strip_item_invalid_inputs():
    with pytest.raises(TypeError):
        strip_item("x", multi_space=0)
    with pytest.raises(TypeError):
        strip_item(42)


# ---------------------------- count_functions --------------------------

def test_count_functions_on_own_module():
    names = count_functions(
        "pycsamt.utils.generic", return_counts=False
    )
    assert "make_ids" in names
    assert "strip_item" in names
    # private and nested helpers excluded by default
    assert all(not n.startswith("_") for n in names)

    n_public = count_functions("pycsamt.utils.generic")
    assert n_public == len(names)
    n_with_private = count_functions(
        "pycsamt.utils.generic", include_private=True
    )
    assert n_with_private >= n_public


def test_count_functions_import_error():
    with pytest.raises(ImportError):
        count_functions("pycsamt.does_not_exist")


# ---------------------------- get_valid_kwargs -------------------------

def test_get_valid_kwargs_filters_and_warns():
    def f(a, b=0, *, c=1):
        return a, b, c

    with pytest.warns(UserWarning):
        valid = get_valid_kwargs(f, {"a": 1, "x": 9, "c": 2})
    assert valid == {"a": 1, "c": 2}


def test_get_valid_kwargs_var_keyword_accepts_all():
    def g(**kw):
        return kw

    assert get_valid_kwargs(g, {"anything": 1}) == {"anything": 1}


def test_get_valid_kwargs_instance_and_class():
    class C:
        def __init__(self, alpha=1, beta=2):
            pass

    assert get_valid_kwargs(C, {"alpha": 5, "nope": 0}.copy()) \
        .keys() == {"alpha"}


# ------------------------------ error_policy ---------------------------

def test_error_policy_passthrough_and_auto():
    assert error_policy("warn") == "warn"
    assert error_policy("raise") == "raise"
    assert error_policy(None, policy="auto", base="warn") == "warn"
    assert error_policy(None, policy=None, base="ignore") == "ignore"


def test_error_policy_strict_and_invalid():
    with pytest.raises(ValueError):
        error_policy(None, policy="strict")
    with pytest.raises(ValueError):
        error_policy("explode")
    with pytest.raises(ValueError):
        error_policy("warn", policy="bogus")
    with pytest.raises(ValueError):
        error_policy(None, policy=None, base="bogus")
    with pytest.raises(KeyError):
        error_policy("explode", exception=KeyError)
