# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""PyInstaller entry script for pyCSAMT Format Studio.

A separate, tiny script (rather than pointing PyInstaller at
``pycsamt/app/converter/__main__.py`` directly) so the frozen build's
entry point is decoupled from the package layout and stays a stable
target for ``pycsamt_converter.spec``.
"""

from pycsamt.app.converter.__main__ import main

if __name__ == "__main__":
    main()
