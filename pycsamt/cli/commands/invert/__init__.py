# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.invert
===========================

Click command group for all inversion-related sub-commands.

Importing this package registers every sub-command on the ``invert``
group.  External code only needs::

    from pycsamt.cli.commands.invert import invert

Sub-command modules
-------------------
_base       Root Click group + shared helpers (solver detection, rich table).
build       ``pycsamt invert build``   — build input files from EDI collection.
run         ``pycsamt invert run``     — launch the inversion binary.
status      ``pycsamt invert status``  — inspect a working directory.
results     ``pycsamt invert results`` — summarise completed results.
plot        ``pycsamt invert plot``    — visualise results (sub-group).
"""

# Importing each sub-module triggers the @invert.command / @invert.group
# decorators, which register the commands on the group defined in _base.
from . import (
    build,  # noqa: F401
    plot,  # noqa: F401
    results,  # noqa: F401
    run,  # noqa: F401
    status,  # noqa: F401
)
from ._base import (
    invert,  # noqa: F401  (re-exported as package public API)
)

__all__ = ["invert"]
