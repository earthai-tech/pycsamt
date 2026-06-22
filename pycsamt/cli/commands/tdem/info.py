# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt tdem info — inspect a TEMAVG survey folder."""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import format_option, no_color_option, verbose_option

from ._base import tdem, _get_survey


@tdem.command("info")
@click.argument(
    "survey_dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="SURVEY_DIR",
)
@click.option(
    "--pattern",
    default="*.AVG",
    show_default=True,
    metavar="GLOB",
    help="Glob pattern used to locate TEMAVG processed-average files.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def info(
    ctx: click.Context,
    survey_dir: Path,
    pattern: str,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Inspect a TEMAVG survey folder and report its contents.

    SURVEY_DIR must be a directory that contains one or more Zonge
    TEMAVG processed-average files (``*.AVG`` by default).

    \b
    Examples:
      pycsamt tdem info data/TEMAVG/JIANGSU
      pycsamt tdem info data/TEMAVG/JIANGSU --format json
      pycsamt tdem info data/TEMAVG/JIANGSU --pattern "*.avg" -v
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    survey = _get_survey(survey_dir, pattern=pattern, verbose=verbose)

    coord = survey.coordinates
    coord_info: dict = {}
    if coord is not None:
        coord_info = {
            "n_points": len(coord.records) if hasattr(coord, "records") else 0,
            "has_elevation": (
                any(getattr(r, "elevation", None) is not None
                    for r in coord.records)
                if hasattr(coord, "records") else False
            ),
        }

    avg_stems = sorted(survey.avg_files)
    z_stems   = sorted(survey.z_files)
    log_stems = sorted(survey.log_files)

    summary = {
        "root": str(survey.root),
        "n_avg_files":  len(avg_stems),
        "n_z_files":    len(z_stems),
        "n_log_files":  len(log_stems),
        "avg_stems":    avg_stems,
        "z_stems":      z_stems,
        "log_stems":    log_stems,
        "has_coordinates": coord is not None,
        "coordinates": coord_info,
    }

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2))
        return

    click.echo(f"Survey root : {survey.root}")
    click.echo(f"AVG files   : {len(avg_stems)}")
    if avg_stems:
        for s in avg_stems:
            avg = survey.avg_files[s]
            n_recs = len(avg.records) if hasattr(avg, "records") else "?"
            click.echo(f"  {s}  ({n_recs} records)")

    click.echo(f"Z files     : {len(z_stems)}")
    if z_stems and verbose >= 1:
        for s in z_stems:
            click.echo(f"  {s}")

    click.echo(f"LOG files   : {len(log_stems)}")
    if log_stems and verbose >= 1:
        for s in log_stems:
            click.echo(f"  {s}")

    if coord is not None:
        n = coord_info.get("n_points", 0)
        has_elev = coord_info.get("has_elevation", False)
        elev_tag = "  (with elevation)" if has_elev else ""
        click.echo(f"Coordinates : {n} points{elev_tag}")
    else:
        click.echo("Coordinates : none found")
