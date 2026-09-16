# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.converter.__main__.

Mirrors ``pycsamt/app/tests/test_desktop_main.py``: PySide6 is imported
lazily inside ``main()``, and ``main()`` instantiates its own
``QApplication`` -- calling it against the module's session-shared
``qapp`` would try to construct a second QApplication in one process, an
unsupported/undefined state. Fake PySide6 modules and a fake
``ConverterMainWindow`` are substituted into ``sys.modules`` instead, so
every real line of ``main()`` runs (icon lookup, ``app.exec()``,
``sys.exit()``) without ever touching real Qt.
"""

from __future__ import annotations

import sys
import types

import pytest

from pycsamt.app.converter import __main__ as main_mod


def test_main_missing_pyside6_exits_with_message(monkeypatch, capsys):
    monkeypatch.setitem(sys.modules, "PySide6.QtWidgets", None)
    monkeypatch.delitem(sys.modules, "PySide6", raising=False)

    with pytest.raises(SystemExit) as exc:
        main_mod.main()

    assert exc.value.code == 1
    err = capsys.readouterr().err
    assert "PySide6 could not be imported for pyCSAMT Format Studio." in err
    assert "pip install 'pycsamt[app]'" in err


def _install_fake_pyside6(monkeypatch):
    app_calls = types.SimpleNamespace(
        argv=None,
        app_name=None,
        display_name=None,
        app_version=None,
        org_name=None,
        icon=None,
        exec_called=False,
    )

    class _FakeQApplication:
        def __init__(self, argv):
            app_calls.argv = argv

        def setApplicationName(self, name):
            app_calls.app_name = name

        def setApplicationDisplayName(self, name):
            app_calls.display_name = name

        def setApplicationVersion(self, version):
            app_calls.app_version = version

        def setOrganizationName(self, name):
            app_calls.org_name = name

        def setWindowIcon(self, icon):
            app_calls.icon = icon

        def exec(self):
            app_calls.exec_called = True
            return 0

    class _FakeQIcon:
        def __init__(self, path):
            self.path = path

    pyside6 = types.ModuleType("PySide6")
    qtwidgets = types.ModuleType("PySide6.QtWidgets")
    qtwidgets.QApplication = _FakeQApplication
    qtgui = types.ModuleType("PySide6.QtGui")
    qtgui.QIcon = _FakeQIcon

    monkeypatch.setitem(sys.modules, "PySide6", pyside6)
    monkeypatch.setitem(sys.modules, "PySide6.QtWidgets", qtwidgets)
    monkeypatch.setitem(sys.modules, "PySide6.QtGui", qtgui)

    return app_calls


def _install_fake_main_window(monkeypatch):
    shown = []

    class _FakeConverterMainWindow:
        def show(self):
            shown.append(True)

    fake_module = types.ModuleType("pycsamt.app.converter.main_window")
    fake_module.ConverterMainWindow = _FakeConverterMainWindow
    monkeypatch.setitem(sys.modules, "pycsamt.app.converter.main_window", fake_module)
    return shown


def test_main_launches_app_successfully(monkeypatch):
    app_calls = _install_fake_pyside6(monkeypatch)
    shown = _install_fake_main_window(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["pycsamt-converter"])

    with pytest.raises(SystemExit) as exc:
        main_mod.main()

    assert exc.value.code == 0
    assert app_calls.argv == ["pycsamt-converter"]
    assert app_calls.app_name == "pycsamt-converter"
    assert app_calls.display_name == "pyCSAMT Format Studio"
    assert app_calls.app_version == "2.0"
    assert app_calls.org_name == "earthai-tech"
    assert app_calls.exec_called is True
    assert shown == [True]
    # The real icon ships under resources/icons/pycsamt.ico.
    assert app_calls.icon is not None


def test_main_skips_icon_when_missing(monkeypatch):
    app_calls = _install_fake_pyside6(monkeypatch)
    _install_fake_main_window(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["pycsamt-converter"])
    monkeypatch.setattr(main_mod.Path, "exists", lambda self: False)

    with pytest.raises(SystemExit):
        main_mod.main()

    assert app_calls.icon is None


def test_module_run_as_script_invokes_main(monkeypatch):
    """Covers the ``if __name__ == "__main__": main()`` guard itself."""
    import runpy

    _install_fake_pyside6(monkeypatch)
    _install_fake_main_window(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["pycsamt-converter"])

    with pytest.raises(SystemExit) as exc:
        runpy.run_module("pycsamt.app.converter.__main__", run_name="__main__")
    assert exc.value.code == 0
