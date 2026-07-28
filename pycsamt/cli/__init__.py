# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli — command-line interface.

Entry point
-----------
``pyproject.toml`` wires ``pycsamt = "pycsamt.cli:main"`` so that
``pip install .`` installs a ``pycsamt`` executable on the PATH.

Programmatic use
----------------
::

    from pycsamt.cli import main

    main(["info", "station.edi"])
    main(["convert", "survey/", "--output-dir", "edis/"])

Available commands
------------------
info        Inspect EDI file / directory metadata.
convert     Convert AVG / J / TEM files to EDI.
"""

from ._base import main

__all__ = ["main"]
