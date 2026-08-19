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

    main(["info", "station.xml"])
    main(["convert", "station.edi", "station.xml"])
    main(["convert", "survey/", "--output-dir", "edis/"])

Available commands
------------------
info        Inspect EDI or EMTF XML transfer-function metadata.
convert     Convert EDI <-> EMTF XML and legacy AVG/J inputs to EDI.
"""

from ._base import main

__all__ = ["main"]
