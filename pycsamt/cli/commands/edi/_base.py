# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.edi._base
================================

Root Click group and shared helpers for ``pycsamt edi``.

The ``edi`` command group is the CLI face of ``pycsamt.seg`` — the
SEG-standard EDI file I/O layer.

Shared helpers
--------------
_get_edi(path, verbose)         → EDIFile  (single file)
_get_collection(path, verbose)  → EDICollection (directory or file list)
"""

from __future__ import annotations

from pathlib import Path

import click


def _get_edi(path: Path, verbose: int = 0):
    """Load a single EDI file and return an :class:`EDIFile`."""
    from pycsamt.seg.edi import EDIFile  # noqa: PLC0415

    edi = EDIFile(path, verbose=verbose)
    edi.read()
    return edi


def _get_collection(path: Path, verbose: int = 0):
    """Load an EDI directory (or file) and return an :class:`EDICollection`."""
    from pycsamt.seg.collection import EDICollection  # noqa: PLC0415

    return EDICollection.from_sources(path, verbose=verbose)


@click.group("edi")
@click.pass_context
def edi(ctx: click.Context) -> None:
    """EDI file I/O — inspect, validate, rotate, filter, and export.

    Operates directly on :class:`~pycsamt.seg.edi.EDIFile` /
    :class:`~pycsamt.seg.collection.EDICollection` objects — the raw
    SEG-standard EDI I/O layer, below the ``Sites`` abstraction used
    by ``pycsamt site``.

    \b
    Sub-commands:
      info       Metadata summary for one EDI file or a whole directory.
      validate   Structural validation — report pass/fail per file.
      stations   Station coordinate table (lat, lon, elev, path).
      profile    Survey profile geometry (bearing, step, distances).
      rotate     Rotate Z and Tipper tensors, write output EDIs.
      select     Filter a collection by station / tipper, export subset.

    \b
    Examples:
      pycsamt edi info    data/AMT/WILLY_DATA/L18PLT/
      pycsamt edi validate data/AMT/ --deep
      pycsamt edi stations data/AMT/WILLY_DATA/L18PLT/ --format csv
      pycsamt edi profile  data/AMT/WILLY_DATA/L18PLT/
      pycsamt edi rotate   data/3edis/ --angle 30 --output-dir rotated/
      pycsamt edi select   data/AMT/WILLY_DATA/L18PLT/ \\
                           --stations 18-023A,18-024A --output-dir subset/
    """
    ctx.ensure_object(dict)
