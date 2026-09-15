# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Unit tests for :mod:`pycsamt.utils.generic`.

``ensure_package`` shells out to pip; its own subprocess.check_call is
mocked below so no real pip invocation ever runs inside the test
suite.
"""

from __future__ import annotations

import subprocess

import numpy as np
import pytest

from pycsamt.utils import generic as generic_mod
from pycsamt.utils.generic import (
    count_functions,
    ensure_package,
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
        "S0",
        "S1",
        "S2",
        "S3",
        "S4",
    ]


def test_make_ids_numpy_count_and_negative():
    assert make_ids(np.int64(2), prefix="n") == ["n0", "n1"]
    with pytest.raises(ValueError):
        make_ids(-1)


# ----------------------------- ensure_package ---------------------------


def test_ensure_package_install_success(monkeypatch):
    calls = {}

    def fake_check_call(cmd, **kw):
        calls["cmd"] = cmd
        calls["kw"] = kw
        return 0

    monkeypatch.setattr(subprocess, "check_call", fake_check_call)
    ok = ensure_package("tqdm", verbose=2, silent=False, capture_output=True)
    assert ok is True
    assert calls["cmd"][:4] == [
        generic_mod.sys.executable, "-m", "pip", "install",
    ]
    assert "--upgrade" in calls["cmd"]


def test_ensure_package_uninstall_and_silent(monkeypatch):
    calls = {}

    def fake_check_call(cmd, **kw):
        calls["cmd"] = cmd
        calls["kw"] = kw
        return 0

    monkeypatch.setattr(subprocess, "check_call", fake_check_call)
    ok = ensure_package(
        "tqdm", install=False, silent=True, extra_pip_args=["--quiet"],
    )
    assert ok is True
    assert "-y" in calls["cmd"]
    assert "--quiet" in calls["cmd"]
    assert calls["kw"]["stdout"] is subprocess.DEVNULL


def test_ensure_package_failure_is_reported_not_raised(monkeypatch):
    def fake_check_call(cmd, **kw):
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(subprocess, "check_call", fake_check_call)
    ok = ensure_package("bogus-package", verbose=1)
    assert ok is False


def test_ensure_package_failure_quiet_when_not_verbose(monkeypatch):
    def fake_check_call(cmd, **kw):
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(subprocess, "check_call", fake_check_call)
    ok = ensure_package("bogus-package", verbose=0)
    assert ok is False


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
    names = count_functions("pycsamt.utils.generic", return_counts=False)
    assert "make_ids" in names
    assert "strip_item" in names
    # private and nested helpers excluded by default
    assert all(not n.startswith("_") for n in names)

    n_public = count_functions("pycsamt.utils.generic")
    assert n_public == len(names)
    n_with_private = count_functions("pycsamt.utils.generic", include_private=True)
    assert n_with_private >= n_public


def test_count_functions_import_error():
    with pytest.raises(ImportError):
        count_functions("pycsamt.does_not_exist")


def test_count_functions_source_unavailable_raises_value_error():
    # Built-in modules have no retrievable Python source.
    with pytest.raises(ValueError, match="Cannot read source"):
        count_functions("sys")


def test_count_functions_include_class_filters_private_classes(tmp_path, monkeypatch):
    module_file = tmp_path / "_count_functions_fixture.py"
    module_file.write_text(
        "def pub_func():\n"
        "    pass\n\n"
        "def _priv_func():\n"
        "    pass\n\n"
        "class PubClass:\n"
        "    pass\n\n"
        "class _PrivClass:\n"
        "    pass\n",
        encoding="utf-8",
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    names = count_functions(
        "_count_functions_fixture", include_class=True, return_counts=False,
    )
    assert names == ["PubClass", "pub_func"]

    count = count_functions("_count_functions_fixture", include_class=True)
    assert count == 2


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

    assert get_valid_kwargs(C, {"alpha": 5, "nope": 0}.copy()).keys() == {"alpha"}


def test_get_valid_kwargs_no_invalid_keys_does_not_warn(recwarn):
    def f(a, b=0):
        return a, b

    valid = get_valid_kwargs(f, {"a": 1, "b": 2})
    assert valid == {"a": 1, "b": 2}
    assert len(recwarn) == 0


def test_get_valid_kwargs_non_callable_instance_resolves_via_class():
    class Plain:
        def __init__(self, alpha=1):
            pass

    instance = Plain()
    assert get_valid_kwargs(instance, {"alpha": 9}) == {"alpha": 9}


def test_get_valid_kwargs_falls_back_to_dunder_call_when_direct_signature_fails(
    monkeypatch,
):
    class Obj:
        def __call__(self, alpha=1, beta=2):
            pass

    obj = Obj()
    real_signature = generic_mod.inspect.signature

    def fake_signature(target):
        if target is obj:
            raise ValueError("no direct signature")
        return real_signature(target)

    monkeypatch.setattr(generic_mod.inspect, "signature", fake_signature)
    result = get_valid_kwargs(obj, {"alpha": 3, "nope": 0})
    assert result == {"alpha": 3}


def test_get_valid_kwargs_dunder_call_also_unresolvable(monkeypatch):
    class Obj:
        def __call__(self, alpha=1):
            pass

    obj = Obj()

    def fake_signature(_target):
        raise TypeError("no signature at all")

    monkeypatch.setattr(generic_mod.inspect, "signature", fake_signature)
    with pytest.warns(UserWarning, match="Unable to retrieve signature"):
        result = get_valid_kwargs(obj, {"alpha": 3})
    assert result == {}


def test_get_valid_kwargs_unresolvable_signature_warns_and_returns_empty(
    monkeypatch,
):
    def _always_raise(_target):
        raise ValueError("no signature")

    monkeypatch.setattr(generic_mod.inspect, "signature", _always_raise)
    with pytest.warns(UserWarning, match="Unable to retrieve signature"):
        result = get_valid_kwargs(object(), {"a": 1})
    assert result == {}


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
