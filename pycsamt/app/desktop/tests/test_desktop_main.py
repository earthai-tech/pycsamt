# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.app.desktop.__main__.

PySide6 is imported lazily inside ``main()`` specifically so this module
stays importable without Qt installed --
``test_main_missing_pyside6_exits_with_message`` exercises that branch and
needs no real PySide6 at all. Every other test here patches
``QApplication``/``QSplashScreen``/``QIcon``/``QPixmap`` directly onto the
real, already-imported ``PySide6.QtWidgets``/``QtGui`` modules (not
whole-module replacement -- PySide6 lazily loads many attributes such as
``QDesktopServices`` via module ``__getattr__`` rather than populating
``__dict__`` up front, so a from-scratch fake module would silently lose
them and break ``LicenseExpiredDialog``/``LicensePage``'s own imports when
the license-gate branch runs). These tests therefore need a real PySide6
install, same as the rest of the desktop test suite.

The license/trial gate (Phase 11) is faked at the ``licensing`` module
level in every test here, never exercised against a real
``OfflineLicenseManager`` -- that would write real trial-tracking state
into this machine's actual QSettings store, an untested-for side effect a
test run should never have.
"""

from __future__ import annotations

import sys
import types

import pytest

from pycsamt.app.desktop import __main__ as desktop_main
from pycsamt.app.desktop.licensing.interfaces import LicenseStatus


def _install_fake_pyside6(
    monkeypatch, pixmap_is_null=False, pixmap_width=200, screen_width=1920
):
    """Inject a minimal fake PySide6 tree so `main()` runs end-to-end
    without a real Qt installation."""
    app_calls = types.SimpleNamespace(
        argv=None,
        app_name=None,
        app_version=None,
        org_name=None,
        org_domain=None,
        icon=None,
        exec_called=False,
        process_events_count=0,
        splash_created_with=None,
        splash_shown=0,
        splash_hidden=0,
        splash_messages=[],
        splash_finished_with=None,
        splash_scaled_to_width=None,
    )

    class _FakeQFont:
        def pointSize(self):
            return 10

        def setPointSize(self, size):
            pass

    class _FakeRect:
        def __init__(self, width):
            self._width = width

        def width(self):
            return self._width

    class _FakeScreen:
        def availableGeometry(self):
            return _FakeRect(screen_width)

    class _FakeQApplication:
        def __init__(self, argv):
            app_calls.argv = argv

        def primaryScreen(self):
            return _FakeScreen()

        def setApplicationName(self, name):
            app_calls.app_name = name

        def setApplicationVersion(self, version):
            app_calls.app_version = version

        def setOrganizationName(self, name):
            app_calls.org_name = name

        def setOrganizationDomain(self, domain):
            app_calls.org_domain = domain

        def setWindowIcon(self, icon):
            app_calls.icon = icon

        def font(self):
            return _FakeQFont()

        def setFont(self, font):
            pass

        def processEvents(self):
            app_calls.process_events_count += 1

        def exec(self):
            app_calls.exec_called = True
            return 0

    class _FakeQIcon:
        def __init__(self, path):
            self.path = path

    class _FakeQPixmap:
        def __init__(self, path):
            self.path = path
            self._width = pixmap_width

        def isNull(self):
            return pixmap_is_null

        def width(self):
            return self._width

        def scaledToWidth(self, width, mode=None):
            scaled = _FakeQPixmap(self.path)
            scaled._width = width
            app_calls.splash_scaled_to_width = width
            return scaled

    class _FakeQSplashScreen:
        def __init__(self, pixmap):
            app_calls.splash_created_with = pixmap

        def show(self):
            app_calls.splash_shown += 1

        def hide(self):
            app_calls.splash_hidden += 1

        def showMessage(self, text, alignment=None, color=None):
            app_calls.splash_messages.append(text)

        def finish(self, widget):
            app_calls.splash_finished_with = widget

    # Built as copy-with-overrides of the real modules, not from scratch:
    # LicenseExpiredDialog/LicensePage (imported when the license-gate
    # branch runs) need real QDialog/QLabel/QLineEdit/QDesktopServices/etc.
    # PySide6 lazily loads many attributes via module __getattr__ rather
    # than populating __dict__ up front, so a from-scratch (or
    # dict-copied) fake module silently loses those -- patch just the
    # handful of names main() itself touches directly on the real,
    # already-imported modules instead of replacing the modules.
    import PySide6.QtGui as _real_qtgui
    import PySide6.QtWidgets as _real_qtwidgets

    monkeypatch.setattr(_real_qtwidgets, "QApplication", _FakeQApplication)
    monkeypatch.setattr(_real_qtwidgets, "QSplashScreen", _FakeQSplashScreen)
    monkeypatch.setattr(_real_qtgui, "QIcon", _FakeQIcon)
    monkeypatch.setattr(_real_qtgui, "QPixmap", _FakeQPixmap)

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


def _install_fake_license_manager(monkeypatch, status=LicenseStatus.TRIAL_ACTIVE):
    """Fake ``licensing.get_default_manager()`` -- never touches real
    QSettings/the registry."""

    class _FakeManager:
        def status(self):
            return status

    manager = _FakeManager()
    monkeypatch.setattr(
        "pycsamt.app.desktop.licensing.get_default_manager", lambda: manager
    )
    return manager


def _install_fake_gate_dialog(monkeypatch, accepted=True):
    """Fake ``LicenseExpiredDialog`` -- records the manager it was built
    with and returns a configurable accept/reject result from ``exec()``."""
    calls = []

    class _DialogCode:
        Rejected = 0
        Accepted = 1

    class _FakeGateDialog:
        DialogCode = _DialogCode

        def __init__(self, manager):
            calls.append(manager)

        def exec(self):
            return _DialogCode.Accepted if accepted else _DialogCode.Rejected

    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.license_expired_dialog.LicenseExpiredDialog",
        _FakeGateDialog,
    )
    return calls


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


def test_main_launches_app_successfully(monkeypatch):
    from pycsamt.app.desktop import branding

    app_calls = _install_fake_pyside6(monkeypatch)
    shown = _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_ACTIVE)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit) as exc:
        desktop_main.main()

    assert exc.value.code == 0
    assert app_calls.argv == ["pycsamt-desktop"]
    assert app_calls.app_name == branding.APP_NAME
    assert app_calls.app_version == branding.get_version()
    assert app_calls.org_name == branding.ORG_NAME
    assert app_calls.org_domain == branding.ORG_DOMAIN
    assert app_calls.exec_called is True
    assert shown == [True]
    # The real icon file ships in resources/icons/pycsamt.logo.ico, so the
    # exists() branch is exercised and setWindowIcon is called.
    assert app_calls.icon is not None


def test_main_skips_icon_when_missing(monkeypatch):
    import pathlib

    app_calls = _install_fake_pyside6(monkeypatch)
    _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_ACTIVE)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])
    monkeypatch.setattr(pathlib.Path, "exists", lambda self: False)

    with pytest.raises(SystemExit):
        desktop_main.main()

    assert app_calls.icon is None
    # The splash image "doesn't exist" either under this blanket patch --
    # confirms the splash is correctly optional, not a hard requirement.
    assert app_calls.splash_created_with is None


def test_module_run_as_script_invokes_main(monkeypatch):
    """Covers the ``if __name__ == "__main__": main()`` guard itself,
    which only executes when the module runs as a script rather than
    being imported (as every other test in this file does)."""
    import runpy

    _install_fake_pyside6(monkeypatch)
    _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_ACTIVE)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit) as exc:
        runpy.run_module("pycsamt.app.desktop.__main__", run_name="__main__")
    assert exc.value.code == 0


# ── Splash screen (Phase 12) ────────────────────────────────────────────────


def test_splash_shown_before_main_window_when_image_exists(monkeypatch):
    """The real splash PNG ships in resources/ -- confirms the happy path
    without needing to fake branding.SPLASH_IMAGE."""
    app_calls = _install_fake_pyside6(monkeypatch)
    _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_ACTIVE)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit):
        desktop_main.main()

    assert app_calls.splash_created_with is not None
    assert app_calls.splash_shown >= 1
    assert "Preparing workspace" in "".join(app_calls.splash_messages)
    assert app_calls.splash_finished_with is not None
    assert app_calls.process_events_count >= 1


def test_splash_skipped_when_pixmap_is_null(monkeypatch):
    """A present-but-corrupt/unreadable image must not crash startup."""
    app_calls = _install_fake_pyside6(monkeypatch, pixmap_is_null=True)
    _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_ACTIVE)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit) as exc:
        desktop_main.main()

    assert exc.value.code == 0
    assert app_calls.splash_shown == 0
    assert app_calls.splash_messages == []


def test_splash_skipped_when_image_path_missing(monkeypatch):
    app_calls = _install_fake_pyside6(monkeypatch)
    _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_ACTIVE)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])
    from pycsamt.app.desktop import branding

    monkeypatch.setattr(
        type(branding.SPLASH_IMAGE), "exists", lambda self: False, raising=False
    )

    with pytest.raises(SystemExit):
        desktop_main.main()

    assert app_calls.splash_shown == 0


def test_oversized_splash_pixmap_is_scaled_down(monkeypatch):
    """Regression test: the real splash PNG ships at 1672x941px (a
    documentation banner image, not a splash-sized asset). Rendered at
    native size, QSplashScreen covers most/all of the screen instead of
    looking like a normal small splash -- caught by the user on a real
    build. Confirms main() scales an oversized pixmap down before
    constructing QSplashScreen."""
    app_calls = _install_fake_pyside6(monkeypatch, pixmap_width=1672, screen_width=1920)
    _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_ACTIVE)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit):
        desktop_main.main()

    assert app_calls.splash_scaled_to_width is not None
    assert app_calls.splash_scaled_to_width <= 480
    assert app_calls.splash_created_with.width() == app_calls.splash_scaled_to_width


def test_splash_scaling_never_exceeds_half_the_available_screen(monkeypatch):
    """A small/narrow screen (e.g. a netbook) must cap the splash width
    below the usual 480px default, not just below the pixmap's own size."""
    app_calls = _install_fake_pyside6(monkeypatch, pixmap_width=1672, screen_width=600)
    _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_ACTIVE)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit):
        desktop_main.main()

    assert app_calls.splash_scaled_to_width == 300  # 600 * 0.5


def test_undersized_splash_pixmap_is_not_scaled(monkeypatch):
    """A pixmap already smaller than the cap must be used as-is -- never
    scaled *up*, which would blur a small icon-sized image."""
    app_calls = _install_fake_pyside6(monkeypatch, pixmap_width=200, screen_width=1920)
    _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_ACTIVE)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit):
        desktop_main.main()

    assert app_calls.splash_scaled_to_width is None
    assert app_calls.splash_created_with.width() == 200


# ── License / trial gate (Phase 11 wiring) ──────────────────────────────────


def test_active_trial_skips_the_gate_dialog_entirely(monkeypatch):
    _install_fake_pyside6(monkeypatch)
    _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_ACTIVE)
    gate_calls = _install_fake_gate_dialog(monkeypatch, accepted=True)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit) as exc:
        desktop_main.main()

    assert exc.value.code == 0
    assert gate_calls == []  # dialog never constructed


def test_licensed_status_skips_the_gate_dialog_entirely(monkeypatch):
    _install_fake_pyside6(monkeypatch)
    _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.LICENSED)
    gate_calls = _install_fake_gate_dialog(monkeypatch, accepted=True)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit) as exc:
        desktop_main.main()

    assert exc.value.code == 0
    assert gate_calls == []


@pytest.mark.parametrize(
    "status", [LicenseStatus.TRIAL_EXPIRED, LicenseStatus.INVALID, LicenseStatus.UNKNOWN]
)
def test_non_active_status_shows_the_gate_dialog(monkeypatch, status):
    app_calls = _install_fake_pyside6(monkeypatch)
    shown = _install_fake_main_window(monkeypatch)
    manager = _install_fake_license_manager(monkeypatch, status)
    gate_calls = _install_fake_gate_dialog(monkeypatch, accepted=True)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit) as exc:
        desktop_main.main()

    assert exc.value.code == 0
    assert gate_calls == [manager]
    # Accepted -> MainWindow still gets built and shown.
    assert shown == [True]
    # Splash is hidden while the gate is up, then re-shown afterward.
    assert app_calls.splash_hidden >= 1


def test_gate_rejected_quits_without_ever_building_main_window(monkeypatch):
    _install_fake_pyside6(monkeypatch)
    shown = _install_fake_main_window(monkeypatch)
    _install_fake_license_manager(monkeypatch, LicenseStatus.TRIAL_EXPIRED)
    gate_calls = _install_fake_gate_dialog(monkeypatch, accepted=False)
    monkeypatch.setattr(sys, "argv", ["pycsamt-desktop"])

    with pytest.raises(SystemExit) as exc:
        desktop_main.main()

    assert exc.value.code == 0
    assert len(gate_calls) == 1
    assert shown == []  # MainWindow never constructed/shown
