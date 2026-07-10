# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt map stations - list or export station coordinates."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    fresh_option,
    no_color_option,
    survey_option,
    verbose_option,
)
from ._base import (
    _format_rows,
    _get_sites,
    _station_rows,
    map,
)


@map.command("stations")
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_DIR",
    required=False,
    default=None,
)
@survey_option
@fresh_option
@click.option(
    "--drop-missing",
    is_flag=True,
    default=False,
    help="Drop stations without finite latitude and longitude.",
)
@click.option(
    "--top",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Print/export only the first N station rows.",
)
@click.option(
    "--output",
    "-o",
    "output_path",
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    default=None,
    metavar="FILE",
    help="Write the rendered table to FILE instead of stdout.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def stations(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    drop_missing: bool,
    top: int | None,
    output_path: Path | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """List station coordinates for a survey.

    EDI_DIR is optional when an active survey has been set with
    ``pycsamt survey set``.  Priority: EDI_DIR > --survey > active context.

    \b
    Examples:
      pycsamt map stations survey/
      pycsamt map stations --drop-missing
      pycsamt map stations --format json
      pycsamt map stations --format csv --output station_coords.csv
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    sites = _get_sites(edi_source, survey_path, fresh, verbose)
    rows = _station_rows(sites, drop_missing=drop_missing)
    if top is not None:
        rows = rows[:top]

    rendered = _format_rows(rows, output_format)

    if output_path is not None:
        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(rendered + "\n", encoding="utf-8")
        except OSError as exc:
            click.echo(f"Error writing {output_path}: {exc}", err=True)
            sys.exit(1)
        click.echo(f"Wrote {len(rows)} station row(s) to {output_path}")
        return

    click.echo(rendered)
