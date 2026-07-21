from __future__ import annotations

import sys
import types

import pytest

from pycsamt.utils._dependency import (
    INSTALL_MAPPING,
    VERSIONS,
    get_version,
    import_optional_dependency,
)

# ------------------------------- get_version ----------------------------


def test_get_version_from_dunder_version():
    m = types.ModuleType("some_mod")
    m.__version__ = "1.2.3"
    assert get_version(m) == "1.2.3"


def test_get_version_from_upper_dunder_version():
    m = types.ModuleType("xlrd")
    m.__VERSION__ = "2.0.1"
    assert get_version(m) == "2.0.1"


def test_get_version_brotli_and_snappy_return_empty():
    brotli_mod = types.ModuleType("brotli")
    snappy_mod = types.ModuleType("snappy")
    assert get_version(brotli_mod) == ""
    assert get_version(snappy_mod) == ""


def test_get_version_unknown_module_raises_import_error():
    m = types.ModuleType("totally_unversioned_mod")
    with pytest.raises(ImportError, match="Can't determine version"):
        get_version(m)


def test_get_version_psycopg2_strips_suffix():
    m = types.ModuleType("psycopg2")
    m.__version__ = "2.9.3 (dt dec pq3 ext lo64)"
    assert get_version(m) == "2.9.3"


# ------------------------- import_optional_dependency --------------------


def test_import_optional_dependency_missing_raises_by_default():
    with pytest.raises(ImportError, match="Missing optional dependency"):
        import_optional_dependency("no_such_pkg_xyz")


def test_import_optional_dependency_missing_ignore_returns_none():
    assert import_optional_dependency("no_such_pkg_xyz", errors="ignore") is None


def test_import_optional_dependency_missing_warn_returns_none():
    assert import_optional_dependency("no_such_pkg_xyz", errors="warn") is None


def test_import_optional_dependency_missing_custom_exception():
    class MyError(Exception):
        pass

    with pytest.raises(MyError):
        import_optional_dependency("no_such_pkg_xyz", exception=MyError)


def test_import_optional_dependency_noncallable_exception_falls_back():
    with pytest.raises(ImportError):
        import_optional_dependency("no_such_pkg_xyz", exception="not_callable")


def test_import_optional_dependency_install_mapping_used_in_message():
    # bs4 maps to beautifulsoup4 in INSTALL_MAPPING and is unlikely installed
    assert INSTALL_MAPPING["bs4"] == "beautifulsoup4"
    try:
        import bs4  # noqa: F401
    except ImportError:
        with pytest.raises(ImportError, match="beautifulsoup4"):
            import_optional_dependency("bs4")
    else:
        pytest.skip("bs4 is installed in this environment")


def test_import_optional_dependency_version_too_old_raises(monkeypatch):
    mod = types.ModuleType("fake_old_pkg")
    mod.__version__ = "0.1"
    monkeypatch.setitem(sys.modules, "fake_old_pkg", mod)
    with pytest.raises(ImportError, match="requires version"):
        import_optional_dependency("fake_old_pkg", min_version="1.0")


def test_import_optional_dependency_version_too_old_warns(monkeypatch):
    mod = types.ModuleType("fake_old_pkg2")
    mod.__version__ = "0.1"
    monkeypatch.setitem(sys.modules, "fake_old_pkg2", mod)
    with pytest.warns(UserWarning, match="requires version"):
        result = import_optional_dependency(
            "fake_old_pkg2", min_version="1.0", errors="warn"
        )
    assert result is None


def test_import_optional_dependency_version_too_old_ignore_returns_module(
    monkeypatch,
):
    mod = types.ModuleType("fake_old_pkg3")
    mod.__version__ = "0.1"
    monkeypatch.setitem(sys.modules, "fake_old_pkg3", mod)
    result = import_optional_dependency(
        "fake_old_pkg3", min_version="1.0", errors="ignore"
    )
    assert result is mod


def test_import_optional_dependency_version_ok_returns_module(monkeypatch):
    mod = types.ModuleType("fake_new_pkg")
    mod.__version__ = "5.0"
    monkeypatch.setitem(sys.modules, "fake_new_pkg", mod)
    result = import_optional_dependency("fake_new_pkg", min_version="1.0")
    assert result is mod


def test_import_optional_dependency_no_min_version_registered(monkeypatch):
    # name not present in VERSIONS -> minimum_version stays None, no check
    mod = types.ModuleType("totally_unregistered_pkg")
    monkeypatch.setitem(sys.modules, "totally_unregistered_pkg", mod)
    assert "totally_unregistered_pkg" not in VERSIONS
    result = import_optional_dependency("totally_unregistered_pkg")
    assert result is mod


def test_import_optional_dependency_submodule_uses_parent_for_version(monkeypatch):
    result = import_optional_dependency("os.path")
    assert result is sys.modules.get("ntpath") or result is sys.modules.get("posixpath")


def test_import_optional_dependency_missing_version_attr_with_min_version(
    monkeypatch,
):
    mod = types.ModuleType("fake_noversion_pkg")
    monkeypatch.setitem(sys.modules, "fake_noversion_pkg", mod)
    with pytest.raises(ImportError, match="Can't determine version"):
        import_optional_dependency("fake_noversion_pkg", min_version="1.0")
