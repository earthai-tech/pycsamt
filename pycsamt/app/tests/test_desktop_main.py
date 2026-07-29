# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.desktop.__main__.

PySide6 is imported lazily inside ``main()`` specifically so this module
stays importable without Qt installed. These tests exercise both branches
by substituting fake modules into ``sys.modules`` rather than requiring a
real PySide6 install, so they run in any environment.
"""

from __future__ import annotations

import sys
import types

import pytest

from pycsamt.app.desktop import __main__ as desktop_main


def test_main_missing_pyside6_exits_with_message(monkeypatch, capsys):
    # Force `from PySide6.QtWidgets import QApplication` to raise ImportError
    # regardless of whether PySide6 is actually installed in this env.
    monkeypatch.setitem(sys.modules, "PySide6.QtWidgets", None)
    monkeypatch.delitem(sys.modules, "PySide6", raising=False)

    with pytest.raises(SystemExit) as exc:
        desktop_main.main()

    assert exc.value.code == 1
    err = capsys.readouterr().err
    assert "PySide6 is required for the desktop app." in err
    assert "pip install 'pycsamt[app]'" in err


def _install_fake_pyside6(monkeypatch):
    """Inject a minimal fake PySide6 tree so `main()` runs end-to-end
    without a real Qt installation."""
    app_calls = types.SimpleNamespace(
        argv=None,
        app_name=None,
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

    class _FakeMainWindow:
        def show(self):
            shown.append(True)

    fake_module = types.ModuleType("pycsamt.app.desktop.main_window")
    fake_module.MainWindow = _FakeMainWindow
    monkeypatch.setitem(sys.modules, "pycsamt.app.desktop.main_window", fake_module)
    return shown


def test_main_launches_app_successfully(monkeypatch):
    app_calls = _install_fake_pyside6(monkeypatch)
    shown = _install_fake_main_window(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit) as exc:
        desktop_main.main()

    assert exc.value.code == 0
    assert app_calls.argv == ["pycsamt-desktop"]
    assert app_calls.app_name == "pycsamt"
    assert app_calls.app_version == "2.0"
    assert app_calls.org_name == "earthai-tech"
    assert app_calls.exec_called is True
    assert shown == [True]
    # The real icon file ships in resources/icons/pycsamt.ico, so the
    # exists() branch is exercised and setWindowIcon is called.
    assert app_calls.icon is not None


def test_main_skips_icon_when_missing(monkeypatch):
    app_calls = _install_fake_pyside6(monkeypatch)
    _install_fake_main_window(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])
    monkeypatch.setattr(desktop_main.Path, "exists", lambda self: False)

    with pytest.raises(SystemExit):
        desktop_main.main()

    assert app_calls.icon is None


def test_module_run_as_script_invokes_main(monkeypatch):
    """Covers the ``if __name__ == "__main__": main()`` guard itself,
    which only executes when the module runs as a script rather than
    being imported (as every other test in this file does)."""
    import runpy

    _install_fake_pyside6(monkeypatch)
    _install_fake_main_window(monkeypatch)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit) as exc:
        runpy.run_module("pycsamt.app.desktop.__main__", run_name="__main__")
    assert exc.value.code == 0
