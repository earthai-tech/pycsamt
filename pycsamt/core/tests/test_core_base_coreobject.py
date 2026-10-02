"""Coverage tests for pycsamt.core.base.CoreObject.

CoreObject backs every core mixin/dataclass (TFBundle, MTBase, BundleMixin)
via its __repr__/__str__/summary/as_dict machinery, but had no direct
tests: this file exercises the field-selection, value-summarization and
dict-coercion logic directly.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest

from pycsamt.core.base import CoreObject


class Plain(CoreObject):
    def __init__(self, name="S001", data=None, _hidden=1):
        self.name = name
        self.data = data if data is not None else [1, 2, 3]
        self._hidden = _hidden


@dataclass
class DC(CoreObject):
    a: int = 1
    b: str = "x"


class NoDict(CoreObject):
    __slots__ = ()


# --------------------------------------------------------------------- #
# __repr__ / __str__ / summary
# --------------------------------------------------------------------- #


def test_repr_shows_public_fields_and_hides_private():
    p = Plain()
    r = repr(p)
    assert r.startswith("Plain(")
    assert "name='S001'" in r
    assert "_hidden" not in r


def test_repr_dataclass_uses_declared_fields():
    d = DC(a=5, b="hello")
    r = repr(d)
    assert r == "DC(a=5, b='hello')"


def test_repr_skips_attribute_that_raises_on_getattr():
    class Weird(CoreObject):
        def __init__(self):
            self.ok = 1

        @property
        def boom(self):
            raise RuntimeError("nope")

    w = Weird()
    # boom is a property, not in __dict__, so _repr_fields (which reads
    # __dict__ keys) never touches it; verify normal fields still show.
    assert "ok=1" in repr(w)


def test_repr_exclude_hook():
    class Excluding(Plain):
        __repr_exclude__ = {"data"}

    e = Excluding()
    r = repr(e)
    assert "name=" in r
    assert "data=" not in r


def test_str_uses_summary():
    p = Plain()
    assert str(p) == p.summary()


def test_summary_truncates_after_max_fields():
    class Many(CoreObject):
        def __init__(self):
            for i in range(10):
                setattr(self, f"f{i}", i)

    m = Many()
    s = m.summary(max_fields=3)
    assert s.count("=") == 3
    assert s.endswith(", ...)")


def test_summary_skips_field_raising_on_access():
    class Flaky(CoreObject):
        def __init__(self):
            self.a = 1
            self.__dict__["b"] = object()

        def _repr_fields(self):
            # force one lookup to raise by yielding a name not truly
            # backed by a real attribute path
            return ["a", "missing_field"]

    f = Flaky()
    s = f.summary()
    assert "a=1" in s
    assert "missing_field" not in s


def test_no_dict_no_dataclass_repr_fields_empty():
    n = NoDict()
    assert list(n._repr_fields()) == []
    assert repr(n) == "NoDict()"


# --------------------------------------------------------------------- #
# _short
# --------------------------------------------------------------------- #


class TestShort:
    def test_primitives(self):
        assert CoreObject._short(None) == "None"
        assert CoreObject._short(True) == "True"
        assert CoreObject._short(3) == "3"
        assert CoreObject._short(3.5) == "3.5"
        assert CoreObject._short(1 + 2j) == "(1+2j)"

    def test_short_string_truncated(self):
        s = "x" * 40
        out = CoreObject._short(s)
        assert out.startswith("'")
        assert "…" in out

    def test_short_string_untruncated(self):
        assert CoreObject._short("hi") == "'hi'"

    def test_short_ndarray(self):
        arr = np.zeros((2, 3), dtype=float)
        out = CoreObject._short(arr)
        assert "ndarray(shape=(2, 3)" in out

    def test_short_mapping(self):
        d = {"a": 1, "b": 2, "c": 3, "d": 4}
        out = CoreObject._short(d)
        assert out.startswith("dict(len=4")
        assert "..." in out

    def test_short_mapping_small(self):
        d = {"a": 1}
        out = CoreObject._short(d)
        assert out == "dict(len=1, keys=[a])"

    def test_short_bytes(self):
        assert CoreObject._short(b"abcd") == "bytes(len=4)"

    def test_short_list_with_sample_and_truncation(self):
        out = CoreObject._short([1, 2, 3, 4, 5])
        assert out.startswith("list([")
        assert "..." in out

    def test_short_empty_list(self):
        assert CoreObject._short([]) == "list([])"

    def test_short_tuple_and_set(self):
        assert CoreObject._short((1, 2)).startswith("tuple([")
        assert CoreObject._short({1, 2}).startswith("set([")

    def test_short_array_like_with_shape(self):
        class ArrLike:
            shape = (4, 4)
            dtype = "float32"

        out = CoreObject._short(ArrLike())
        assert "array_like(shape=(4, 4)" in out

    def test_short_fallback_repr(self):
        class Custom:
            def __repr__(self):
                return "<Custom!>"

        assert CoreObject._short(Custom()) == "<Custom!>"


# --------------------------------------------------------------------- #
# as_dict / _to_dict / _coerce
# --------------------------------------------------------------------- #


def test_as_dict_dataclass():
    d = DC(a=7, b="y")
    out = d.as_dict()
    assert out == {"a": 7, "b": "y"}


def test_as_dict_generic_object_public_only():
    p = Plain(name="S1", data=[1, 2])
    out = p.as_dict()
    assert out["name"] == "S1"
    assert out["data"] == [1, 2]
    assert "_hidden" not in out


def test_as_dict_generic_object_include_private():
    p = Plain(name="S1")
    out = p.as_dict(public_only=False)
    assert "_hidden" in out


def test_as_dict_max_depth_zero_returns_marker():
    p = Plain()
    out = p.as_dict(max_depth=-1)
    assert out == {"_": "max_depth"}


def test_as_dict_nested_object_one_level():
    class Inner(CoreObject):
        def __init__(self):
            self.x = 1

    class Outer(CoreObject):
        def __init__(self):
            self.inner = Inner()

    out = Outer().as_dict(max_depth=2)
    assert out["inner"] == {"x": 1}


def test_coerce_ndarray_summarized():
    class Holder(CoreObject):
        def __init__(self):
            self.arr = np.arange(6).reshape(2, 3)

    out = Holder().as_dict()
    assert out["arr"]["type"] == "ndarray"
    assert out["arr"]["shape"] == (2, 3)


def test_coerce_mapping_limits_to_32_items():
    class Holder(CoreObject):
        def __init__(self):
            self.d = {i: i for i in range(40)}

    out = Holder().as_dict()
    assert len(out["d"]) == 32


def test_coerce_sequence_limits_to_32_items():
    class Holder(CoreObject):
        def __init__(self):
            self.seq = list(range(40))

    out = Holder().as_dict()
    assert len(out["seq"]) == 32


def test_object_without_dict_to_dict_returns_empty():
    class Slotted(CoreObject):
        __slots__ = ("x",)

        def __init__(self):
            self.x = 1

    out = CoreObject._to_dict(Slotted(), public_only=True, depth=1)
    assert out == {}
