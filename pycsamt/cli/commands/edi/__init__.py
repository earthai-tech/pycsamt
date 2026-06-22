# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.edi
=========================

Click command group for EDI file I/O operations.

This is the CLI face of ``pycsamt.seg`` — the SEG-standard EDI I/O layer.
It operates directly on :class:`~pycsamt.seg.edi.EDIFile` and
:class:`~pycsamt.seg.collection.EDICollection`, below the ``Sites``
abstraction used by ``pycsamt site``.

Sub-command modules
-------------------
_base       Root Click group + ``_get_edi`` / ``_get_collection`` helpers.
info        ``pycsamt edi info``      — metadata summary (file or directory).
validate    ``pycsamt edi validate``  — structural EDI validation.
stations    ``pycsamt edi stations``  — station coordinate table.
profile     ``pycsamt edi profile``   — survey profile geometry.
rotate      ``pycsamt edi rotate``    — rotate Z / Tipper, write EDIs.
select      ``pycsamt edi select``    — filter collection, export subset.
"""

from ._base import edi  # noqa: F401

from . import info      # noqa: F401  (registers @edi.command("info"))
from . import validate  # noqa: F401  (registers @edi.command("validate"))
from . import stations  # noqa: F401  (registers @edi.command("stations"))
from . import profile   # noqa: F401  (registers @edi.command("profile"))
from . import rotate    # noqa: F401  (registers @edi.command("rotate"))
from . import select    # noqa: F401  (registers @edi.command("select"))

__all__ = ["edi"]
