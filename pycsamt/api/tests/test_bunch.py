from __future__ import annotations

import pytest

from pycsamt.api.bunch import Bunch, FlexDict


# ─────────────────────────────────────────────────────────────────────────
# Bunch
# ─────────────────────────────────────────────────────────────────────────


def test_bunch_attribute_access():
    b = Bunch(a=1, b=2)
    assert b.a == 1
    assert b["b"] == 2


def test_bunch_getattr_missing_raises_attribute_error():
    b = Bunch(a=1)
    with pytest.raises(AttributeError):
        b.missing


def test_bunch_setattr_sets_item():
    b = Bunch()
    b.x = 42
    assert b["x"] == 42


def test_bunch_setstate_is_a_noop():
    b = Bunch(a=1)
    b.__setstate__({"whatever": "ignored"})
    assert dict(b) == {"a": 1}


def test_bunch_dir_includes_keys():
    b = Bunch(a=1, b=2)
    assert "a" in dir(b)
    assert "b" in dir(b)


def test_bunch_repr_lists_keys():
    b = Bunch(a=1, b=2)
    text = repr(b)
    assert text.startswith("<Bunch object with keys:")
    assert "a" in text and "b" in text


def test_bunch_str_renders_without_dunder_keys():
    b = Bunch(a=1, __hidden__=2)
    text = str(b)
    assert isinstance(text, str)
    assert "a" in text


def test_bunch_str_empty_bunch():
    b = Bunch()
    text = str(b)
    assert isinstance(text, str)


# ─────────────────────────────────────────────────────────────────────────
# FlexDict
# ─────────────────────────────────────────────────────────────────────────


def test_flexdict_attribute_access():
    fd = FlexDict(a=1, b=2)
    assert fd.a == 1
    assert fd["b"] == 2


def test_flexdict_getattr_missing_raises_attribute_error():
    fd = FlexDict(a=1)
    with pytest.raises(AttributeError):
        fd.missing


def test_flexdict_setattr_sets_item():
    fd = FlexDict()
    fd.x = 42
    assert fd["x"] == 42


def test_flexdict_setattr_strips_special_symbols():
    fd = FlexDict()
    fd.__setattr__("name**suffix", "value")
    assert fd["name"] == "value"
    assert "name**suffix" not in fd


@pytest.mark.parametrize("symbol", ["**", "%%", "&&", "||", "$$"])
def test_flexdict_setattr_strips_each_special_symbol(symbol):
    fd = FlexDict()
    fd.__setattr__(f"base{symbol}rest", "v")
    assert fd["base"] == "v"


def test_flexdict_setstate_updates_and_rebinds_dict():
    fd = FlexDict(a=1)
    fd.__setstate__({"b": 2})
    assert fd["a"] == 1
    assert fd["b"] == 2
    assert fd.__dict__ is fd
    assert "__dict__" not in fd


def test_flexdict_dir_returns_keys_only():
    fd = FlexDict(a=1, b=2)
    assert set(dir(fd)) == {"a", "b"}


def test_flexdict_repr_lists_keys():
    fd = FlexDict(a=1, b=2)
    text = repr(fd)
    assert text.startswith("<FlexDict with keys:")
    assert "a" in text and "b" in text
