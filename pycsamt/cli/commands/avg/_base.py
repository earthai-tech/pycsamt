# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.avg._base
================================

Root Click group and shared helpers for ``pycsamt avg``.

The ``avg`` command group is the CLI face of ``pycsamt.zonge`` — the
Zonge CSAMT/AMT AVG file processing ecosystem.

Shared helpers
--------------
_get_avg(path, verbose)  → AVG object (modern K2 path uses AVG.from_file;
                           legacy K1 path falls back to raw loader)
_load_raw(path)          → (df, meta, kind) without xarray requirement
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click


def _get_avg(path: Path, verbose: int = 0):
    """Load a Zonge AVG file and return an :class:`~pycsamt.zonge.avg.AVG`.

    Tries ``AVG.from_file()`` first.  When the file is legacy (kind-1)
    and *xarray* is not installed the function falls back to the raw
    DataFrame loader and raises :class:`click.UsageError` with a helpful
    message so the user knows what is needed.
    """
    from pycsamt.zonge.avg import AVG  # noqa: PLC0415

    try:
        return AVG.from_file(path, verbose=bool(verbose))
    except ImportError as exc:
        raise click.UsageError(
            f"Loading {path.name!r} requires an optional dependency:\n\n"
            f"  {exc}\n\n"
            "Install it with:  pip install xarray\n\n"
            "Modern (kind-2) AVG files do not need xarray."
        ) from exc


def _load_raw(path: Path):
    """Return *(df, meta, kind)* from an AVG file without xarray.

    Uses ``load_avg()`` and ``classify_avg_format()`` directly.
    Suitable for QC and metadata inspection when xarray is absent.
    """
    from pycsamt.zonge.utils import load_avg, classify_avg_format  # noqa: PLC0415

    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    kind = classify_avg_format(lines)
    df, meta = load_avg(path)
    return df, meta, kind


@click.group("avg")
@click.pass_context
def avg(ctx: click.Context) -> None:
    """Zonge AVG file I/O — inspect, validate QC, correct, and export.

    Operates on Zonge CSAMT/AMT processed-average (``.avg``) files using
    :class:`~pycsamt.zonge.avg.AVG` (modern kind-2) or the raw loader
    (legacy kind-1).  For AVG → EDI conversion use ``pycsamt transform avg``.

    \b
    Sub-commands:
      info       Metadata summary (project, hardware, stations, frequencies).
      validate   QC table — % errors, phase sigma, SNR, skip flags.
      stations   Station positions from an optional .stn coordinate file.
      correct    Apply static-shift or capacitive-coupling correction.
      export     Write processed data to modern, legacy, or netCDF format.

    \b
    Examples:
      pycsamt avg info    data/avg/K2.AVG
      pycsamt avg validate data/avg/K2.AVG --format csv
      pycsamt avg stations data/avg/K2.AVG --stn-file data/avg/K2.stn
      pycsamt avg correct  data/avg/K2.AVG --method tma --ref-freq 2048
      pycsamt avg export   data/avg/K2.AVG --format modern --output-dir out/
    """
    ctx.ensure_object(dict)
