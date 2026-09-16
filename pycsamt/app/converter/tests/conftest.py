# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Shared fixtures for pycsamt.app.converter tests.

This directory is a sibling of ``pycsamt/app/tests`` (not a
subdirectory), so pytest does not pick up ``pycsamt/app/tests/conftest.py``
automatically -- the Qt offscreen setup, shared ``qapp``, and modal-dialog
neutralization needed by GUI tests here are duplicated (trimmed to what
this suite actually needs) rather than imported, since conftest.py files
are not importable modules across sibling packages.
"""

from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

_qt_active = False


@pytest.fixture(scope="session", autouse=True)
def qt_offscreen():
    """Force Qt to use the offscreen (headless) platform for all tests."""
    global _qt_active
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    try:
        import PySide6  # noqa: F401
    except ImportError:
        yield
        return
    _qt_active = True
    yield


@pytest.fixture(scope="session")
def qapp():
    """Single QApplication shared across the whole test session."""
    from PySide6.QtWidgets import QApplication

    global _qt_active
    _qt_active = True

    existing = QApplication.instance()
    if existing is not None:
        yield existing
        return

    app = QApplication(["pytest", "-platform", "offscreen"])
    app.setApplicationName("pycsamt-converter-test")
    app.setOrganizationName("earthai-tech")
    yield app


@pytest.fixture(autouse=True)
def _settle_qt_gc_after_test():
    """Reclaim cyclic Qt garbage right after each test, not mid the next.

    See pycsamt/app/tests/conftest.py for the full rationale (same
    Shiboken teardown-ordering hazard applies to any Qt widget test).
    """
    yield
    if _qt_active:
        from pycsamt.compat.qt import settle_qt_gc

        settle_qt_gc()


@pytest.fixture(autouse=True)
def no_global_restyle(monkeypatch):
    """Make ``QApplication.setStyleSheet`` a no-op (see app/tests/conftest.py)."""
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        yield
        return

    monkeypatch.setattr(QApplication, "setStyleSheet", lambda self, *_a, **_k: None)
    yield


@pytest.fixture(autouse=True)
def no_modal_dialogs(monkeypatch):
    """Neutralize modal QMessageBox dialogs (see app/tests/conftest.py)."""
    try:
        from PySide6.QtWidgets import QMessageBox
    except ImportError:
        yield
        return

    yes = QMessageBox.StandardButton.Yes
    ok = QMessageBox.StandardButton.Ok
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: yes))
    monkeypatch.setattr(QMessageBox, "information", staticmethod(lambda *a, **k: ok))
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda *a, **k: ok))
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a, **k: ok))
    yield


@pytest.fixture()
def isolated_qsettings(tmp_path, monkeypatch):
    """Redirect ConverterSettings' QSettings backend to a throwaway ini file.

    Without this, ``load_settings``/``save_settings`` hit the real
    per-user registry/config store (``earthai-tech`` / ``pycsamt-converter``)
    and tests would leak state into the developer's actual environment.
    """
    from PySide6.QtCore import QSettings

    from pycsamt.app.converter import settings as settings_mod

    ini_path = tmp_path / "converter_settings.ini"

    def _fake_qsettings():
        return QSettings(str(ini_path), QSettings.Format.IniFormat)

    monkeypatch.setattr(settings_mod, "_qsettings", _fake_qsettings)
    return ini_path
