# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Entry point: python -m pycsamt.app.desktop"""

import sys
from pathlib import Path


def main() -> None:
    # PySide6 import deferred so the module is importable without Qt installed.
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError as exc:
        print(
            "PySide6 is required for the desktop app.\n"
            "Install it with:  pip install 'pycsamt[app]'",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc

    from PySide6.QtGui import QIcon
    from pycsamt.app.desktop.main_window import MainWindow

    app = QApplication(sys.argv)
    app.setApplicationName("pycsamt")
    app.setApplicationVersion("2.0")
    app.setOrganizationName("earthai-tech")

    # Application icon — multi-resolution ICO
    _ico = Path(__file__).parent / "resources" / "icons" / "pycsamt.ico"
    if _ico.exists():
        app.setWindowIcon(QIcon(str(_ico)))

    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
