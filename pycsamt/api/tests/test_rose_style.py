from __future__ import annotations

import pytest

from pycsamt.api._rose_style import RoseStyle, resolve_rose_style


def test_rose_style_copy_applies_overrides_without_mutating_original():
    rs = RoseStyle()
    copied = rs.copy(compass_labels="degrees", show_secondary=False)
    assert copied.compass_labels == "degrees"
    assert copied.show_secondary is False
    assert rs.compass_labels == "NESW"
    assert rs.show_secondary is True


def test_rose_style_copy_rejects_unknown_attribute():
    rs = RoseStyle()
    with pytest.raises(ValueError, match="has no attribute"):
        rs.copy(bogus_attr=1)


def test_rose_style_repr_lists_public_fields():
    rs = RoseStyle()
    text = repr(rs)
    assert text.startswith("RoseStyle(")
    assert "bar_style=" in text


def test_resolve_rose_style_none_returns_default_preset():
    rs = resolve_rose_style(None)
    assert rs.cmap == "YlOrRd"


def test_resolve_rose_style_named_presets():
    minimal = resolve_rose_style("minimal")
    assert minimal.show_mean is False
    publication = resolve_rose_style("PUBLICATION")
    assert publication.compass_labels == "degrees"
    alias = resolve_rose_style("pycsamt_rose")
    assert alias.cmap == "YlOrRd"


def test_resolve_rose_style_unknown_name_raises():
    with pytest.raises(ValueError, match="Unknown rose style"):
        resolve_rose_style("not-a-real-preset")


def test_resolve_rose_style_instance_passthrough():
    custom = RoseStyle(cmap="Blues")
    assert resolve_rose_style(custom) is custom


def test_resolve_rose_style_instance_with_overrides_returns_copy():
    custom = RoseStyle(cmap="Blues")
    resolved = resolve_rose_style(custom, compass_labels="degrees")
    assert resolved is not custom
    assert resolved.compass_labels == "degrees"
    assert resolved.cmap == "Blues"


def test_resolve_rose_style_rejects_bad_type():
    with pytest.raises(TypeError, match="must be a str, RoseStyle, or None"):
        resolve_rose_style(42)
