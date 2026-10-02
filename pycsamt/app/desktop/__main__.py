# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Entry point: python -m pycsamt.app.desktop"""

import os

# Several desktop features (AI Inversion, AI-processing agents --
# pycsamt.ai.* / pycsamt.agents.*) import torch/tensorflow/sklearn, and
# other desktop code paths (Occam2D forward modeling, curve fitting) use
# real scipy.optimize -- both in the same long-lived process. On
# Windows/conda that combination can abort the whole interpreter ("OMP:
# Error #15: Initializing libiomp5md.dll, but found ... already
# initialized") -- duplicate OpenMP runtime registration, not a bug in
# either library. Same root cause already documented and guarded against
# in every test conftest that imports both (see MEMORY.md); the desktop
# app itself never had this guard. Must be set before either is imported.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import sys
import warnings


def main() -> None:
    # PySide6 import deferred so the module is importable without Qt installed.
    try:
        from PySide6.QtCore import Qt, qInstallMessageHandler
        from PySide6.QtWidgets import QApplication
    except ImportError as exc:
        print(
            "PySide6 is required for the desktop app.\n"
            "Install it with:  pip install 'pycsamt[app]'",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc

    from PySide6.QtGui import QIcon, QPixmap
    from PySide6.QtWidgets import QSplashScreen

    from pycsamt.app.desktop import branding

    # Qt can emit one harmless invalid-point-size message while discovering
    # the Windows system font, before QApplication lets us correct that font.
    previous_qt_handler = None

    def _qt_message_handler(msg_type, context, message):
        if "QFont::setPointSize: Point size <= 0" in message:
            return
        if previous_qt_handler is not None:
            previous_qt_handler(msg_type, context, message)
        else:
            print(message, file=sys.stderr)

    previous_qt_handler = qInstallMessageHandler(_qt_message_handler)
    app = QApplication(sys.argv)
    # Some Windows themes expose only a pixel-sized default font.  Ensure
    # widgets never try to derive an invalid point size (-1) from it.
    app_font = app.font()
    if app_font.pointSize() <= 0:
        app_font.setPointSize(10)
        app.setFont(app_font)

    # Matplotlib may report cosmetic tight-layout limitations asynchronously
    # from the Qt event loop. Embedded canvases use explicit margins instead.
    warnings.filterwarnings(
        "ignore",
        message=r"Tight layout not applied\..*",
        category=UserWarning,
    )
    warnings.filterwarnings(
        "ignore",
        message=r"This figure includes Axes that are not compatible with tight_layout.*",
        category=UserWarning,
    )
    # Same for constrained layout: a canvas drawn while its panel is still
    # a few pixels wide (hidden window, collapsed splitter, first paint)
    # cannot fit its axes; it lays out again at its real size.
    warnings.filterwarnings(
        "ignore",
        message=r"constrained_layout not applied because axes sizes "
        r"collapsed to zero.*",
        category=UserWarning,
    )
    app.setApplicationName(branding.APP_NAME)
    app.setApplicationVersion(branding.get_version())
    app.setOrganizationName(branding.ORG_NAME)
    app.setOrganizationDomain(branding.ORG_DOMAIN)

    # Application icon — multi-resolution ICO
    if branding.LOGO_ICO.exists():
        app.setWindowIcon(QIcon(str(branding.LOGO_ICO)))

    # ── Splash screen (Phase 12) ────────────────────────────────────────
    # Shown before any heavy import (MainWindow pulls in matplotlib/pandas/
    # every panel controller) so something appears on screen immediately
    # instead of a multi-second blank delay after launch.
    splash = None
    if branding.SPLASH_IMAGE.exists():
        pixmap = QPixmap(str(branding.SPLASH_IMAGE))
        if not pixmap.isNull():
            # The source PNG is a 1672x941 documentation banner image, not
            # a splash-sized asset -- QSplashScreen renders a pixmap at its
            # native size with no auto-scaling, so left as-is this covers
            # most or all of the screen instead of looking like a normal
            # small splash. Cap it to a sane splash width, scaled down
            # (never up) and never wider than half the available screen.
            max_width = 480
            screen = app.primaryScreen()
            if screen is not None:
                max_width = min(max_width, int(screen.availableGeometry().width() * 0.5))
            if pixmap.width() > max_width:
                pixmap = pixmap.scaledToWidth(
                    max_width, Qt.TransformationMode.SmoothTransformation
                )
            splash = QSplashScreen(pixmap)
            splash.show()
            app.processEvents()

    # ── License / trial gate (Phase 11) ─────────────────────────────────
    # Hard block, per the product decision in the modernization plan's
    # Section 4: MainWindow is never constructed without an active trial
    # or a valid license.
    from pycsamt.app.desktop.licensing import LicenseStatus, get_default_manager

    manager = get_default_manager()
    if manager.status() not in (LicenseStatus.TRIAL_ACTIVE, LicenseStatus.LICENSED):
        if splash is not None:
            splash.hide()

        from pycsamt.app.desktop.dialogs.license_expired_dialog import (
            LicenseExpiredDialog,
        )

        gate = LicenseExpiredDialog(manager)
        if gate.exec() != gate.DialogCode.Accepted:
            sys.exit(0)

        if splash is not None:
            splash.show()
            app.processEvents()

    # ── Preparing workspace (cosmetic step, not a dependency installer —
    # a PyInstaller-frozen build ships every dependency already bundled;
    # see the modernization plan's Guiding Principle 3) ─────────────────
    if splash is not None:
        splash.showMessage(
            "Preparing workspace…",
            alignment=Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignHCenter,
            color=Qt.GlobalColor.white,
        )
        app.processEvents()

    # Heavy import, deferred until after the splash/license gate above so
    # the splash is what the user sees first, not a blank window.
    from pycsamt.app.desktop.main_window import MainWindow

    window = MainWindow()
    window.show()
    if splash is not None:
        splash.finish(window)
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
