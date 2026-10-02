# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt airborne — root Click group and shared helpers.

Every sub-module imports ``airborne`` from here and registers its
command via ``@airborne.command(...)``.  The ``__init__.py`` imports
each sub-module to trigger the registrations.

Shared utility
--------------
_get_asites(source, recursive, verbose) -> AirborneSites
    Coerce a path (single EMTF-XML file or directory) to an
    :class:`~pycsamt.airborne.site.AirborneSites` collection.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import click

if TYPE_CHECKING:
    from pycsamt.airborne.site import AirborneSites


def _get_asites(
    source: Path,
    recursive: bool,
    verbose: int,
) -> AirborneSites:
    """Coerce *source* (EMTF-XML file or directory) to ``AirborneSites``."""
    from pycsamt.airborne.site import ensure_asites  # noqa: PLC0415

    if verbose >= 1:
        click.echo(f"Loading airborne sites from {source}", err=True)
    return ensure_asites(source, recursive=recursive)


@click.group("airborne")
@click.pass_context
def airborne(ctx: click.Context) -> None:
    """Inspect and diagnose airborne EM datasets (ZTEM, MobileMT, AFMAG).

    Operates on already-decoded EMTF-XML files or directories of them
    (a single site file, or a whole flight survey) -- pyCSAMT has no
    raw-vendor-file reader for airborne systems, only generic adapters.

    \b
    Typical workflow:
      pycsamt airborne info      survey_xml/
      pycsamt airborne diagnose  survey_xml/
      pycsamt airborne diagnose  survey_xml/ --format csv --output table.csv
    """
    ctx.ensure_object(dict)
