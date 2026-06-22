# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.jones
===========================

Click command group for Jones J-format file operations.

Operates on :class:`~pycsamt.jones.j.JFile` and
:class:`~pycsamt.jones.collection.JCollection` — the native A.G. Jones
J-format I/O layer.

For J → EDI conversion see ``pycsamt transform j``.

Sub-command modules
-------------------
_base      Root Click group + ``_get_jfile`` / ``_get_collection`` helpers.
info       ``pycsamt jones info``      — metadata summary (file or directory).
validate   ``pycsamt jones validate``  — structural J-file validation.
stations   ``pycsamt jones stations``  — station coordinate and metadata table.
blocks     ``pycsamt jones blocks``    — data block inspection (kind/comp/rows).
select     ``pycsamt jones select``    — filter collection, export subset.
"""

from ._base import jones  # noqa: F401

from . import info      # noqa: F401  (registers @jones.command("info"))
from . import validate  # noqa: F401  (registers @jones.command("validate"))
from . import stations  # noqa: F401  (registers @jones.command("stations"))
from . import blocks    # noqa: F401  (registers @jones.command("blocks"))
from . import select    # noqa: F401  (registers @jones.command("select"))

__all__ = ["jones"]
