from __future__ import annotations

import re

import pytest

from pycsamt.utils.handlers import columns_manager

# ---------------------------- columns_manager -----------------------------


def test_columns_manager_none_default_empty_as_none_true():
    assert columns_manager(None) is None


def test_columns_manager_none_empty_as_none_false():
    assert columns_manager(None, empty_as_none=False) == []


def test_columns_manager_none_with_default():
    assert columns_manager(None, default=["a", "b"]) == ["a", "b"]


def test_columns_manager_string_default_pattern():
    # default pattern only splits on [@&,;#], whitespace is preserved
    assert columns_manager("col1@col2&col3") == ["col1", "col2", "col3"]


def test_columns_manager_string_custom_separator():
    assert columns_manager("col1,col2,col3", separator=",") == [
        "col1",
        "col2",
        "col3",
    ]


def test_columns_manager_string_custom_regex():
    regex = re.compile(r"\s*\|\s*")
    assert columns_manager("col1 | col2 | col3", regex=regex) == [
        "col1",
        "col2",
        "col3",
    ]


def test_columns_manager_string_custom_pattern_kwarg():
    assert columns_manager("col1#col2#col3", pattern=r"#") == [
        "col1",
        "col2",
        "col3",
    ]


def test_columns_manager_list_and_tuple_passthrough():
    assert columns_manager(["a", "b"]) == ["a", "b"]
    assert columns_manager(("a", "b")) == ["a", "b"]


def test_columns_manager_numeric_scalar_wrapped():
    assert columns_manager(5) == [5]
    assert columns_manager(5.5) == [5.5]


def test_columns_manager_callable_wrapped():
    def f():
        pass

    result = columns_manager(f)
    assert result == [f]


def test_columns_manager_dict_wrap_true():
    d = {"a": 1, "b": 2}
    assert columns_manager(d, wrap_dict=True) == [d]


def test_columns_manager_dict_wrap_false_uses_keys():
    d = {"a": 1, "b": 2}
    assert columns_manager(d, wrap_dict=False) == ["a", "b"]


def test_columns_manager_class_object_wrapped():
    class Foo:
        pass

    assert columns_manager(Foo) == [Foo]


def test_columns_manager_non_iterable_instance_wrapped():
    class Foo:
        pass

    instance = Foo()
    assert columns_manager(instance) == [instance]


def test_columns_manager_to_upper_all_strings():
    assert columns_manager(["a", "b"], to_upper=True) == ["A", "B"]


def test_columns_manager_to_upper_mixed_types_raise():
    with pytest.raises(TypeError):
        columns_manager([1, "b"], to_upper=True, error="raise")


def test_columns_manager_to_upper_mixed_types_warn():
    with pytest.warns(UserWarning, match="skipping 'upper'"):
        result = columns_manager([1, "b"], to_upper=True, error="warn")
    assert result == [1, "b"]


def test_columns_manager_to_upper_mixed_types_ignore():
    result = columns_manager([1, "b"], to_upper=True, error="ignore")
    assert result == [1, "b"]


def test_columns_manager_to_string():
    assert columns_manager([1, 2], to_string=True) == ["1", "2"]


class _BadIterOnce:
    """Iterable whose iteration raises, so list() conversion fails."""

    def __iter__(self):
        def gen():
            raise RuntimeError("boom")
            yield  # pragma: no cover - unreachable, makes this a generator

        return gen()


def test_columns_manager_list_conversion_error_raise():
    with pytest.raises(ValueError, match="Error converting columns to list"):
        columns_manager(_BadIterOnce(), error="raise")


def test_columns_manager_list_conversion_error_warn():
    with pytest.warns(UserWarning, match="Could not convert columns to list"):
        result = columns_manager(_BadIterOnce(), error="warn")
    assert isinstance(result, _BadIterOnce)


def test_columns_manager_list_conversion_error_ignore():
    obj = _BadIterOnce()
    result = columns_manager(obj, error="ignore")
    assert result is obj


def test_columns_manager_generic_iterable_converted_to_list():
    assert columns_manager(range(3)) == [0, 1, 2]


class _BadIterImmediate:
    """Object whose __iter__ call itself raises (not just the generator)."""

    def __iter__(self):
        raise RuntimeError("no iter")


def test_columns_manager_final_fallback_wraps_when_iter_also_fails():
    obj = _BadIterImmediate()
    result = columns_manager(obj, error="ignore")
    assert result == [obj]


class _BadTuple(tuple):
    """Tuple subclass whose __iter__ raises, tripping the list-conversion
    error path while still satisfying isinstance(columns, tuple)."""

    def __iter__(self):
        raise RuntimeError("bad tuple iter")


def test_columns_manager_tuple_subclass_survives_failed_conversion():
    bt = _BadTuple((1, 2, 3))
    result = columns_manager(bt, error="ignore")
    assert result is bt


class _InstanceLevelIter:
    """hasattr(obj, '__iter__') is True (instance attribute) while
    isinstance(obj, Iterable) is False (ABC check looks at the class)."""


def test_columns_manager_hasattr_iter_but_not_abc_iterable_wrapped():
    obj = _InstanceLevelIter()
    obj.__iter__ = lambda: iter([1, 2, 3])
    result = columns_manager(obj)
    assert result == [obj]
