from __future__ import annotations

import threading
from dataclasses import dataclass

import numpy as np
import pytest

from pycsamt.api.property import MetadataMixin, PyCSAMTObject


# ─────────────────────────────────────────────────────────────────────────
# __repr__ / __str__ / summary
# ─────────────────────────────────────────────────────────────────────────


class _RaisingAttr(PyCSAMTObject):
    __repr_fields__ = ("a", "boom")

    def __init__(self):
        self.a = 1

    @property
    def boom(self):
        raise RuntimeError("nope")


def test_repr_skips_attributes_that_raise():
    obj = _RaisingAttr()
    text = repr(obj)
    assert "a=1" in text
    assert "boom" not in text


class _ManyFields(PyCSAMTObject):
    __repr_max_fields__ = 2

    def __init__(self):
        self.a = 1
        self.b = 2
        self.c = 3


def test_repr_truncates_with_ellipsis_when_over_max_fields():
    text = repr(_ManyFields())
    assert "..." in text


class _CustomStr(PyCSAMTObject):
    def summary(self, *, max_fields=None):
        return 123  # not a string


def test_str_falls_back_to_repr_when_summary_not_str():
    obj = _CustomStr()
    assert str(obj) == repr(123)


def test_summary_with_max_fields_restores_original_after():
    obj = _ManyFields()
    original = obj.__repr_max_fields__
    text = obj.summary(max_fields=1)
    assert obj.__repr_max_fields__ == original
    assert isinstance(text, str)


def test_summary_without_max_fields_returns_repr():
    obj = _ManyFields()
    assert obj.summary() == repr(obj)


# ─────────────────────────────────────────────────────────────────────────
# update / clone / validate
# ─────────────────────────────────────────────────────────────────────────


class _Validated(PyCSAMTObject):
    def __init__(self, x=0):
        self.x = x
        self.validated = False

    def validate(self):
        self.validated = True
        if self.x < 0:
            raise ValueError("x must be non-negative")


def test_update_sets_attributes_and_validates():
    obj = _Validated()
    result = obj.update(x=5)
    assert result is obj
    assert obj.x == 5
    assert obj.validated is True


def test_update_propagates_validation_errors():
    obj = _Validated()
    with pytest.raises(ValueError):
        obj.update(x=-1)


def test_clone_returns_shallow_copy_with_overrides():
    obj = _Validated(x=1)
    clone = obj.clone(x=9)
    assert clone is not obj
    assert clone.x == 9
    assert obj.x == 1
    assert clone.validated is True


def test_base_validate_is_a_noop():
    class Plain(PyCSAMTObject):
        pass

    assert Plain().validate() is None


# ─────────────────────────────────────────────────────────────────────────
# _slot_names: single-string __slots__
# ─────────────────────────────────────────────────────────────────────────


class _SingleSlot(PyCSAMTObject):
    __slots__ = "only_slot"

    def __init__(self):
        self.only_slot = "value"


def test_slot_names_handles_bare_string_slots():
    obj = _SingleSlot()
    assert obj.to_dict() == {"only_slot": "value"}


# ─────────────────────────────────────────────────────────────────────────
# _short_value branches
# ─────────────────────────────────────────────────────────────────────────


def test_short_value_truncates_long_strings():
    long_string = "x" * 100
    out = PyCSAMTObject._short_value(long_string, max_string=10)
    assert out.endswith("...'")
    assert len(out) < len(repr(long_string))


def test_short_value_ndarray():
    arr = np.zeros((2, 3), dtype=float)
    out = PyCSAMTObject._short_value(arr)
    assert "ndarray(shape=(2, 3)" in out


def test_short_value_bytes():
    out = PyCSAMTObject._short_value(b"abc")
    assert out == "bytes(len=3)"


def test_short_value_sequence_with_ellipsis_when_long():
    out = PyCSAMTObject._short_value([1, 2, 3, 4, 5])
    assert "..." in out
    assert "len=5" in out


def test_short_value_mapping_with_few_keys_has_no_ellipsis():
    out = PyCSAMTObject._short_value({"a": 1, "b": 2})
    assert "..." not in out
    assert "len=2" in out


def test_short_value_mapping_with_many_keys_has_ellipsis():
    out = PyCSAMTObject._short_value({"a": 1, "b": 2, "c": 3, "d": 4})
    assert "..." in out
    assert "len=4" in out


def test_short_value_array_like_without_ndarray_type():
    class FakeArray:
        shape = (3, 3)
        dtype = "float64"

    out = PyCSAMTObject._short_value(FakeArray())
    assert "array_like(shape=(3, 3)" in out


def test_short_value_falls_back_to_repr():
    class Plain:
        def __repr__(self):
            return "<plain>"

    assert PyCSAMTObject._short_value(Plain()) == "<plain>"


# ─────────────────────────────────────────────────────────────────────────
# to_dict / _to_dict branches
# ─────────────────────────────────────────────────────────────────────────


def test_to_dict_max_depth_reached():
    assert PyCSAMTObject._to_dict(
        {"a": 1}, public_only=True, depth=-1
    ) == {"_": "max_depth"}


def test_to_dict_ndarray():
    out = PyCSAMTObject._to_dict(
        np.ones(3), public_only=True, depth=1,
    )
    assert out == {"type": "ndarray", "shape": (3,), "dtype": "float64"}


@dataclass
class _PlainDataclass:
    a: int
    b: str


_PlainDataclass.__repr_exclude__ = {"b"}


def test_to_dict_dataclass_uses_asdict_and_excludes_fields():
    obj = _PlainDataclass(a=1, b="secret")
    out = PyCSAMTObject._to_dict(obj, public_only=True, depth=1)
    assert out == {"a": 1}


@dataclass
class _UndeepcopyableDataclass:
    lock: object
    skip_me: int = 0


_UndeepcopyableDataclass.__repr_exclude__ = {"skip_me"}


def test_to_dict_dataclass_falls_back_when_asdict_fails():
    obj = _UndeepcopyableDataclass(lock=threading.Lock(), skip_me=1)
    out = PyCSAMTObject._to_dict(obj, public_only=True, depth=1)
    assert set(out.keys()) == {"lock"}


def test_to_dict_mapping_limits_to_32_and_recurses():
    mapping = {f"k{i}": i for i in range(40)}
    out = PyCSAMTObject._to_dict(mapping, public_only=True, depth=1)
    assert len(out) == 32


def test_to_dict_sequence_limits_to_32_and_recurses():
    out = PyCSAMTObject._to_dict(
        list(range(40)), public_only=True, depth=1,
    )
    assert len(out) == 32
    assert out[0] == 0


def test_to_dict_generic_object_public_only_excludes_private_and_excluded():
    class Obj(PyCSAMTObject):
        __repr_exclude__ = {"hidden"}

        def __init__(self):
            self.visible = 1
            self._private = 2
            self.hidden = 3

    out = PyCSAMTObject._to_dict(Obj(), public_only=True, depth=1)
    assert out == {"visible": 1}


def test_to_dict_generic_object_includes_private_when_not_public_only():
    class Obj(PyCSAMTObject):
        def __init__(self):
            self.visible = 1
            self._private = 2

    out = PyCSAMTObject._to_dict(Obj(), public_only=False, depth=1)
    assert out == {"visible": 1, "_private": 2}


def test_to_dict_generic_object_skips_attributes_that_raise():
    class Obj(PyCSAMTObject):
        def __init__(self):
            self.ok = 1

        @property
        def boom(self):
            raise RuntimeError("nope")

    class Sub(Obj):
        __dict__ = {}  # placeholder so hasattr check works normally

    obj = Obj()
    obj.__dict__["boom"] = None  # ensure "boom" appears as a name to try
    # Overwrite with a descriptor-raising property lookup via getattr
    out = PyCSAMTObject._to_dict(obj, public_only=True, depth=1)
    assert "ok" in out


# ─────────────────────────────────────────────────────────────────────────
# MetadataMixin
# ─────────────────────────────────────────────────────────────────────────


class _WithMetadata(MetadataMixin):
    pass


def test_metadata_mixin_ensure_creates_dict_when_missing():
    obj = _WithMetadata()
    meta = obj.ensure_metadata()
    assert meta == {}
    assert obj.metadata is meta


def test_metadata_mixin_update_and_snapshot():
    obj = _WithMetadata()
    result = obj.update_metadata(a=1, b=2)
    assert result is obj
    snapshot = obj.metadata_dict()
    assert snapshot == {"a": 1, "b": 2}
    snapshot["a"] = 999
    assert obj.metadata["a"] == 1  # snapshot is a shallow copy
