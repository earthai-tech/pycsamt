# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt tdem — root Click group and shared helpers.

Every sub-module imports ``tdem`` from here and registers its commands
via ``@tdem.command(...)`` or ``@tdem.group(...)``.  The ``__init__.py``
imports each sub-module to trigger the registrations.

Shared utility
--------------
_get_survey(path, pattern, verbose) → TEMSurvey
    Parse a TEMAVG survey folder and return a :class:`TEMSurvey`.
"""

from __future__ import annotations

from pathlib import Path

import click


def _get_survey(path: Path, pattern: str = "*.AVG", verbose: int = 0):
    """Parse a TEMAVG survey folder and return a :class:`TEMSurvey`.

    Raises :class:`click.UsageError` when no AVG files are found.
    """
    from pycsamt.tdem.survey import read_temavg_survey  # noqa: PLC0415

    survey = read_temavg_survey(path, pattern=pattern, verbose=verbose)
    if not survey.avg_files:
        raise click.UsageError(
            f"No TEMAVG files matching {pattern!r} found in {path}.\n\n"
            "  Check that the directory contains Zonge TEMAVG .AVG files."
        )
    return survey


@click.group("tdem")
@click.pass_context
def tdem(ctx: click.Context) -> None:
    """Time-domain EM (TEMAVG) loading, conversion, and visualisation.

    \b
    Typical workflow:
      pycsamt tdem info    data/TEMAVG/JIANGSU
      pycsamt tdem convert data/TEMAVG/JIANGSU --output-dir imp/
      pycsamt tdem plot    data/TEMAVG/JIANGSU --kind decay

    \b
    Sub-commands:
      info       Inspect a TEMAVG survey folder.
      convert    Convert TEMAVG → frequency-domain EDI files.
      plot       Visualise decay curves, apparent-ρ sections, or maps.
    """
    ctx.ensure_object(dict)
