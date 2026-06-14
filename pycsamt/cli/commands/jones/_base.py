# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.jones._base
=================================

Root Click group and shared helpers for ``pycsamt jones``.

Shared helpers
--------------
_get_jfile(path, verbose)       → JFile  (single file)
_get_collection(path, verbose)  → JCollection (directory or file list)

The A.G. Jones J-format stores MT data as one text file per station
with a banner, an ``>KEY=VALUE`` info block, and homogeneous data blocks
(R for rho/phase, Z for complex impedance, T/C for GDS transfer functions).
"""

from __future__ import annotations

from pathlib import Path

import click


def _get_jfile(path: Path, verbose: int = 0):
    """Load a single Jones J-file and return a :class:`JFile`."""
    from pycsamt.jones.j import JFile  # noqa: PLC0415

    return JFile.from_file(path, verbose=verbose)


def _get_collection(path: Path, verbose: int = 0):
    """Load a directory (or list) of J-files and return a :class:`JCollection`."""
    from pycsamt.jones.collection import JCollection  # noqa: PLC0415

    return JCollection.from_sources(path, verbose=verbose)


@click.group("jones")
@click.pass_context
def jones(ctx: click.Context) -> None:
    """Jones J-format file inspection, validation, filtering, and export.

    Operates on :class:`~pycsamt.jones.j.JFile` /
    :class:`~pycsamt.jones.collection.JCollection` objects — the native
    A.G. Jones J-format I/O layer.

    For J → EDI conversion see ``pycsamt transform j``.

    \b
    Sub-commands:
      info       Metadata summary for one J-file or a directory.
      validate   Structural validation — pass/fail per file.
      stations   Station table (lat, lon, elev, azimuth, n_freq).
      blocks     Inspect data blocks in a J-file (kind, component, rows).
      select     Filter a collection and export a subset.

    \b
    Examples:
      pycsamt jones info    data/j/
      pycsamt jones validate data/j/ --deep
      pycsamt jones stations data/j/ --format csv
      pycsamt jones blocks   data/j/S01.j
      pycsamt jones select   data/j/ --has-z --output-dir subset/
    """
    ctx.ensure_object(dict)
