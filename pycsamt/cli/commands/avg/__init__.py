# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.avg
=========================

Click command group for Zonge CSAMT/AMT AVG file operations.

This is the CLI face of ``pycsamt.zonge`` — the Zonge processed-average
file ecosystem.  For AVG → EDI conversion use ``pycsamt transform avg``.

Sub-command modules
-------------------
_base       Root Click group + ``_get_avg`` / ``_load_raw`` helpers.
info        ``pycsamt avg info``      — metadata summary.
validate    ``pycsamt avg validate``  — QC table (% errors, phase sigma).
stations    ``pycsamt avg stations``  — station positions and coordinates.
correct     ``pycsamt avg correct``   — static-shift / capacitive-coupling.
export      ``pycsamt avg export``    — write modern, legacy, or netCDF.
"""

from . import (
    correct,  # noqa: F401  (registers @avg.command("correct"))
    export,  # noqa: F401  (registers @avg.command("export"))
    info,  # noqa: F401  (registers @avg.command("info"))
    stations,  # noqa: F401  (registers @avg.command("stations"))
    validate,  # noqa: F401  (registers @avg.command("validate"))
)
from ._base import avg  # noqa: F401

__all__ = ["avg"]
