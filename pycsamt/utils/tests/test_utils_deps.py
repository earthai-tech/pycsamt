from __future__ import annotations

import sys
import types

import pytest

from pycsamt.utils.deps import ensure_pkg

# ------------------------------- ensure_pkg -------------------------------


def test_ensure_pkg_present_package_wraps_function():
    decorator = ensure_pkg("os")

    def func(x):
        """doc"""
        return x + 1

    wrapped = decorator(func)
    assert wrapped(1) == 2
    assert wrapped.__name__ == "func"
    assert wrapped.__doc__ == "doc"


def test_ensure_pkg_missing_raises_by_default():
    with pytest.raises(ImportError):
        ensure_pkg("no_such_pkg_xyz")


def test_ensure_pkg_missing_errors_ignore_does_not_raise():
    decorator = ensure_pkg("no_such_pkg_xyz", errors="ignore")

    def func():
        return 42

    assert decorator(func)() == 42


def test_ensure_pkg_auto_install_success_retries_and_succeeds(monkeypatch):
    calls = {}

    def fake_ensure_package(dist_name, install=True, upgrade=True, verbose=0):
        calls["dist_name"] = dist_name
        calls["verbose"] = verbose
        mod = types.ModuleType("fake_autoinstalled_pkg")
        mod.__version__ = "1.0"
        sys.modules["fake_autoinstalled_pkg"] = mod
        return True

    monkeypatch.setattr("pycsamt.utils.deps.ensure_package", fake_ensure_package)
    try:
        decorator = ensure_pkg("fake_autoinstalled_pkg", auto_install=True, verbose=2)
        assert calls["dist_name"] == "fake_autoinstalled_pkg"

        def func():
            return "ok"

        assert decorator(func)() == "ok"
    finally:
        sys.modules.pop("fake_autoinstalled_pkg", None)


def test_ensure_pkg_auto_install_failure_reraises_original_error(monkeypatch):
    def fake_ensure_package(dist_name, install=True, upgrade=True, verbose=0):
        return False

    monkeypatch.setattr("pycsamt.utils.deps.ensure_package", fake_ensure_package)
    with pytest.raises(ImportError):
        ensure_pkg("no_such_pkg_xyz", auto_install=True)


def test_ensure_pkg_auto_install_still_missing_after_install(monkeypatch):
    # ensure_package reports success but the module still can't be imported
    def fake_ensure_package(dist_name, install=True, upgrade=True, verbose=0):
        return True

    monkeypatch.setattr("pycsamt.utils.deps.ensure_package", fake_ensure_package)
    result_decorator = ensure_pkg("no_such_pkg_xyz", auto_install=True, errors="ignore")

    def func():
        return "still ok"

    assert result_decorator(func)() == "still ok"


def test_ensure_pkg_auto_install_uses_dist_name_override(monkeypatch):
    seen = {}

    def fake_ensure_package(dist_name, install=True, upgrade=True, verbose=0):
        seen["dist_name"] = dist_name
        mod = types.ModuleType("fake_pkg_with_dist_name")
        mod.__version__ = "1.0"
        sys.modules["fake_pkg_with_dist_name"] = mod
        return True

    monkeypatch.setattr("pycsamt.utils.deps.ensure_package", fake_ensure_package)
    try:
        ensure_pkg(
            "fake_pkg_with_dist_name",
            dist_name="some-other-pypi-name",
            auto_install=True,
        )
        assert seen["dist_name"] == "some-other-pypi-name"
    finally:
        sys.modules.pop("fake_pkg_with_dist_name", None)


def test_ensure_pkg_auto_install_min_version_upgrade_path(monkeypatch):
    mod = types.ModuleType("fake_old_version_pkg")
    mod.__version__ = "0.1"
    monkeypatch.setitem(sys.modules, "fake_old_version_pkg", mod)

    def fake_ensure_package(dist_name, install=True, upgrade=True, verbose=0):
        mod.__version__ = "9.0"
        return True

    monkeypatch.setattr("pycsamt.utils.deps.ensure_package", fake_ensure_package)
    decorator = ensure_pkg("fake_old_version_pkg", min_version="1.0", auto_install=True)

    def func():
        return "upgraded"

    assert decorator(func)() == "upgraded"


def test_ensure_pkg_auto_install_verbose_zero_no_print(monkeypatch, capsys):
    def fake_ensure_package(dist_name, install=True, upgrade=True, verbose=0):
        return False

    monkeypatch.setattr("pycsamt.utils.deps.ensure_package", fake_ensure_package)
    with pytest.raises(ImportError):
        ensure_pkg("no_such_pkg_xyz", auto_install=True, verbose=0)
    captured = capsys.readouterr()
    assert captured.err == ""
