"""Backward-compat shim: pycsamt.emtools._rose_style re-exports api._rose_style."""

from pycsamt.api._rose_style import (
    _PRESETS,
    _UNSET,
    RoseStyle,
    resolve_rose_style,
)
from pycsamt.emtools import _rose_style


def test_reexports_match_canonical_module():
    assert _rose_style._PRESETS is _PRESETS
    assert _rose_style._UNSET is _UNSET
    assert _rose_style.RoseStyle is RoseStyle
    assert _rose_style.resolve_rose_style is resolve_rose_style
