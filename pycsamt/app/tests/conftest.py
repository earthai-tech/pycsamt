# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Shared fixtures for pycsamt.app tests.

Qt is initialised once per session with the offscreen platform so
tests can instantiate Qt widgets without a real display.
"""

from __future__ import annotations

import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# Some tests import both torch and scipy in the same process (e.g.
# pycsamt.app.web.callbacks.inversion pulls in torch at module scope, then
# a "Traditional" inversion test calls scipy.optimize.least_squares for
# real). On Windows/conda that combination can abort the whole interpreter
# ("OMP: Error #15: Initializing libiomp5md.dll, but found ... already
# initialized") -- duplicate OpenMP runtime registration, not a bug in
# either library. Must be set before either one is imported.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

_exit_status = 0
_qt_active = False
_IS_WINDOWS = sys.platform == "win32"


def pytest_sessionfinish(session, exitstatus):  # noqa: ARG001
    global _exit_status
    _exit_status = int(exitstatus)


def pytest_unconfigure(config):  # noqa: ARG001
    """
    Hard-kill the process to dodge a PySide6/Shiboken shutdown crash.

    Widgets created across the whole test session build up reference
    cycles that CPython's normal shutdown GC resolves in an order
    Shiboken doesn't expect (``QObject: shared QObject was deleted
    directly``), segfaulting *after* results and coverage have already
    been written. Nothing meaningful happens between here and process
    exit, so leave before interpreter shutdown touches the Qt object
    graph — without importing Qt or querying ``QApplication.instance()``
    here, since that state may already be corrupted.

    ``os._exit`` maps straight to the ``_exit(2)`` syscall on POSIX
    (where CI actually runs): it skips ``Py_Finalize`` entirely, and
    since shared objects loaded via the dynamic linker get no
    equivalent of a DLL-unload notification on ``_exit``, that's
    enough. On Windows, ``ExitProcess`` (which ``os._exit`` still goes
    through) *does* run ``DLL_PROCESS_DETACH`` for every loaded DLL,
    including Qt's native ones, so it crashes the same way there;
    ``TerminateProcess`` is the only call that skips that too, and it
    takes the exit code directly so status reporting stays intact.
    """
    if not _qt_active:
        return
    sys.stdout.flush()
    sys.stderr.flush()
    if _IS_WINDOWS:
        import ctypes

        kernel32 = ctypes.windll.kernel32
        kernel32.TerminateProcess(kernel32.GetCurrentProcess(), _exit_status)
    else:
        os._exit(_exit_status)


# ── Qt / offscreen setup ───────────────────────────────────────────────────


@pytest.fixture(scope="session", autouse=True)
def qt_offscreen():
    """Force Qt to use the offscreen (headless) platform for all tests."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
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
    app.setApplicationName("pycsamt-test")
    yield app


@pytest.fixture(scope="session", autouse=True)
def _stub_mpl_qt_toolbar():
    """
    Replace matplotlib's Qt ``NavigationToolbar2QT`` with a lightweight stub.

    Constructing the real toolbar (``MplCanvas(toolbar=True)`` ->
    ``NavigationToolbar2QT(canvas, self)``) under the offscreen Qt platform
    while coverage's C tracer is active corrupts interpreter memory on the
    Python 3.9 / matplotlib 3.9.x combination and segfaults -- the crash
    surfaces mid-construction, e.g. in ``pathlib`` while the toolbar loads its
    icons (CI: "Python 3.9 / interfaces" shard, test_map_window.py). It is the
    same 3.9-only fragility that the terminal matplotlib 3.9.x line shows
    elsewhere; Python >=3.10 pulls matplotlib >=3.10 and is unaffected.

    The toolbar is a pure matplotlib convenience widget carrying no pycsamt
    logic, and no test asserts on it, so swapping it for a ``QWidget`` that
    satisfies the little the widget touches (``.canvas`` attribute,
    ``sizeHint()``, ``update()``) removes the crash without losing meaningful
    coverage. ``MplCanvas`` resolves the class lazily via ``_qt_mpl_classes``
    (``from matplotlib.backends.backend_qtagg import NavigationToolbar2QT``),
    so patching the name on that module is enough.

    Scoped to Python < 3.10; set ``PYCSAMT_STUB_MPL_TOOLBAR=1`` to force it on
    newer interpreters (used to verify the stub is a valid drop-in locally).
    """
    force = os.environ.get("PYCSAMT_STUB_MPL_TOOLBAR") == "1"
    if not (force or sys.version_info < (3, 10)):
        yield
        return
    try:
        from matplotlib.backends import backend_qtagg
        from PySide6.QtWidgets import QWidget
    except Exception:
        yield
        return

    class _StubNavigationToolbar(QWidget):
        def __init__(self, canvas, parent=None):
            super().__init__(parent)
            self.canvas = canvas

    original = backend_qtagg.NavigationToolbar2QT
    backend_qtagg.NavigationToolbar2QT = _StubNavigationToolbar
    try:
        yield
    finally:
        backend_qtagg.NavigationToolbar2QT = original


@pytest.fixture(autouse=True)
def no_global_restyle(monkeypatch):
    """
    Make ``QApplication.setStyleSheet`` a no-op.

    Applying an application-wide stylesheet forces Qt to repolish
    every live widget (including matplotlib canvases), which hangs
    under the offscreen platform. No test asserts on the stylesheet
    text; theme tests assert on session state and icons instead.
    """
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        yield
        return

    monkeypatch.setattr(
        QApplication, "setStyleSheet", lambda self, *_a, **_k: None
    )
    yield


@pytest.fixture(autouse=True)
def no_modal_dialogs(monkeypatch):
    """
    Neutralize modal QMessageBox dialogs.

    Modal dialogs block forever under the offscreen platform (e.g.
    the quit confirmation in ``MainWindow.closeEvent``), hanging the
    whole suite. Every static box returns its affirmative button.
    """
    try:
        from PySide6.QtWidgets import QMessageBox
    except ImportError:
        yield
        return

    yes = QMessageBox.StandardButton.Yes
    ok = QMessageBox.StandardButton.Ok
    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(lambda *a, **k: yes),
    )
    monkeypatch.setattr(
        QMessageBox,
        "information",
        staticmethod(lambda *a, **k: ok),
    )
    monkeypatch.setattr(
        QMessageBox,
        "warning",
        staticmethod(lambda *a, **k: ok),
    )
    monkeypatch.setattr(
        QMessageBox,
        "critical",
        staticmethod(lambda *a, **k: ok),
    )
    yield


# ── Dash web app ─────────────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def web_app():
    """Single fully-wired Dash app (all callbacks registered) for reuse."""
    from pycsamt.app.web.app import create_app

    return create_app()


# ── EDI fixtures ───────────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def simulated_edi(tmp_path_factory):
    """Return a path to a minimal synthetic EDI file."""
    import numpy as np

    tmp = tmp_path_factory.mktemp("edi")
    edi_path = tmp / "SIM001.edi"

    # Build the minimal >INFO / >HEAD / >FREQ / >ZXXR … blocks
    nfreq = 8
    freqs = np.logspace(2, -1, nfreq)
    z_real = np.ones(nfreq) * 10.0
    np.ones(nfreq) * 8.0

    lines = [
        ">HEAD",
        " DATAID=SIM001",
        " LAT=48:30:0.0",
        " LONG=7:45:0.0",
        " ELEV=200.0",
        ">INFO",
        " MAXINFO=999",
        ">DEFINEMEAS",
        " MAXCHAN=7",
        " MAXRUN=999",
        " MAXMEAS=9999",
        " UNITS=M",
        " REFTYPE=CART",
        " REFLAT=48:30:0.0",
        " REFLONG=7:45:0.0",
        " REFELEV=200.0",
        ">=MTSECT",
        " SECTID=SIM001",
        f" NFREQ={nfreq}",
        " HX=1001.001",
        " HY=1002.001",
        " HZ=1003.001",
        " EX=1004.001",
        " EY=1005.001",
        f">FREQ // {nfreq}",
    ]
    lines.append("  " + "  ".join(f"{f:.6E}" for f in freqs))
    for comp in (
        "ZXXR",
        "ZXXI",
        "ZXYR",
        "ZXYI",
        "ZYXR",
        "ZYXI",
        "ZYYR",
        "ZYYI",
    ):
        lines.append(f">{comp} // {nfreq}")
        lines.append("  " + "  ".join(f"{v:.6E}" for v in z_real))
    lines.append(">END")

    edi_path.write_text("\n".join(lines), encoding="utf-8")
    return edi_path
