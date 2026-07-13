from __future__ import annotations

import types
import warnings

import pytest

from pycsamt.compat.aliases import (
    compat_alias,
    install_compat_aliases,
    make_compat_alias,
)


def _target_a(x, y=1, *, z=0):
    """Simple target for tests."""
    return x + y + z


def _target_b(v):
    return v * 2


def test_make_compat_alias_warn_and_forward():
    alias = make_compat_alias(
        _target_a,
        old_name="old_a",
        new_name="new_a",
        since="2.0.0",
        remove_in="2.1.0",
        extra="Read the v2 guide.",
    )

    # emits FutureWarning and forwards the call
    with pytest.warns(FutureWarning) as rec:
        out = alias(2, 3, z=5)
    assert out == _target_a(2, 3, z=5)

    # message content checks (substring to keep robust)
    msg = str(rec[0].message)
    assert "old_a" in msg
    assert "new_a" in msg
    assert "deprecated since v2.0.0" in msg
    assert "removed in v2.1.0" in msg
    assert "Read the v2 guide." in msg

    # marker attributes and preserved name
    assert getattr(alias, "__is_compat_alias__", False)
    assert getattr(alias, "__compat_target__", None) is _target_a
    assert alias.__name__ == "old_a"


def test_install_compat_aliases_basic_and_skip():
    ns: dict[str, object] = {}

    mapping = {
        "old_a": _target_a,
        "old_b": _target_b,
    }
    extras = {"old_b": "Use keyword 'style' in v2."}

    # first install
    install_compat_aliases(
        mapping,
        g=ns,
        since="2.0.0",
        remove_in="3.0.0",
        extras=extras,
    )

    assert "old_a" in ns and "old_b" in ns

    # calls warn and forward
    with pytest.warns(FutureWarning) as rec_a:
        out_a = ns["old_a"](1, 2, z=3)
    assert out_a == _target_a(1, 2, z=3)
    assert "deprecated since v2.0.0" in str(rec_a[0].message)

    # Be extra-robust against env filters on FutureWarning.
    with warnings.catch_warnings(record=True) as rec_b:
        warnings.simplefilter("always", FutureWarning)
        out_b = ns["old_b"](7)

    assert out_b == _target_b(7)
    assert any(issubclass(w.category, FutureWarning) for w in rec_b)
    assert any("style" in str(w.message) for w in rec_b)

    # prepare an existing alias; second install must skip it
    existing = make_compat_alias(
        _target_a,
        old_name="old_a",
        new_name="new_a",
        since="2.0.0",
        remove_in="3.0.0",
    )
    ns["old_a"] = existing

    install_compat_aliases(
        mapping,
        g=ns,
        since="2.0.0",
        remove_in="3.0.0",
        extras=extras,
    )

    # still the same object (skipped)
    assert ns["old_a"] is existing


def test_compat_alias_injects_and_exports(monkeypatch):
    # create a clean module-like namespace
    mod = types.ModuleType("tmp_mod")
    mod.__dict__["__all__"] = ["new_func"]

    @compat_alias(
        "old_func",
        since="2.0.0",
        remove_in="3.0.0",
        extra="See v2 docs.",
        export=True,
    )
    def new_func(a, b=0):
        return a - b

    # the decorator injects into the defining module, which is
    # this test module's globals, not 'mod'. To validate export
    # behavior safely, re-apply inside 'mod' via exec.
    src = (
        "from pycsamt.compat.aliases import compat_alias\n"
        "__all__ = ['new_func']\n"
        "@compat_alias('old_func', remove_in='3.0.0')\n"
        "def new_func(a, b=0):\n"
        "    return a - b\n"
    )
    exec(src, mod.__dict__, mod.__dict__)

    # alias exists and warns; function itself does not warn
    assert "old_func" in mod.__dict__
    assert "new_func" in mod.__dict__

    with pytest.warns(FutureWarning):
        out = mod.__dict__["old_func"](5, b=2)
    assert out == 3

    # direct call: no warning
    # with pytest.warns(None) as rec_none:
    #     out2 = mod.__dict__["new_func"](5, b=2)
    # Ensure direct call to the new function does not warn.
    with warnings.catch_warnings(record=True) as rec_none:
        warnings.simplefilter("always", FutureWarning)
        out2 = mod.__dict__["new_func"](5, b=2)

    assert out2 == 3
    assert len(rec_none) == 0

    # export appended to __all__
    assert "old_func" in mod.__dict__["__all__"]


if __name__ == "__main__":  # pragma: no-cover
    pytest.main([__file__])
