# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Entry point: python -m pycsamt.app.converter"""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    # PySide6 import deferred so the module is importable without Qt installed.
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError as exc:
        print(
            "PySide6 could not be imported for pyCSAMT Format Studio.\n"
            f"  {type(exc).__name__}: {exc}\n"
            "If PySide6 is not installed, get it with:  "
            "pip install 'pycsamt[app]'",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc

    from PySide6.QtGui import QIcon

    from pycsamt.app.converter.main_window import ConverterMainWindow

    app = QApplication(sys.argv)
    app.setApplicationName("pycsamt-converter")
    app.setApplicationDisplayName("pyCSAMT Format Studio")
    app.setApplicationVersion("2.0")
    app.setOrganizationName("earthai-tech")

    icon_path = Path(__file__).parent / "resources" / "icons" / "pycsamt.ico"
    if icon_path.exists():
        app.setWindowIcon(QIcon(str(icon_path)))

    # Theme (light/dark) stylesheet is applied by ConverterMainWindow itself,
    # from the persisted setting -- same division of responsibility as the
    # full desktop app's __main__.py / MainWindow.
    window = ConverterMainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
