# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt map plot - save a static station-location map."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    fresh_option,
    no_color_option,
    survey_option,
    verbose_option,
)

from ._base import _get_sites, _station_rows, map


@map.command("plot")
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
    "--output",
    "-o",
    "output_path",
    type=click.Path(dir_okay=False, writable=True, path_type=Path),
    default=Path("station_map.png"),
    show_default=True,
    metavar="FILE",
    help="Figure file to write.",
)
@click.option(
    "--title",
    default="Station Map",
    show_default=True,
    metavar="TEXT",
    help="Figure title.",
)
@click.option(
    "--label/--no-label",
    default=True,
    show_default=True,
    help="Annotate station markers with station names.",
)
@click.option(
    "--dpi",
    type=click.IntRange(min=50),
    default=150,
    show_default=True,
    help="Figure resolution in dots per inch.",
)
@click.option(
    "--show",
    is_flag=True,
    default=False,
    help="Display an interactive matplotlib window after saving.",
)
@verbose_option
@no_color_option
@click.pass_context
def plot(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    output_path: Path,
    title: str,
    label: bool,
    dpi: int,
    show: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Plot station longitude/latitude positions as a static figure.

    This first map command is intentionally dependency-light: it uses
    matplotlib and the station coordinates already stored in the EDI
    headers.  Rich basemap and interactive exports can be layered on this
    command group later without changing the survey-resolution API.

    \b
    Examples:
      pycsamt map plot survey/ --output station_map.png
      pycsamt map plot --no-label --dpi 200
      pycsamt map plot --show
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    sites = _get_sites(edi_source, survey_path, fresh, verbose)
    rows = _station_rows(sites, drop_missing=True)
    if not rows:
        raise click.UsageError(
            "No stations with finite latitude and longitude were found."
        )

    try:
        from pycsamt.map import (  # noqa: PLC0415
            StationMapOptions,
            build_station_map,
            ensure_map_data,
        )
        data = ensure_map_data(sites, verbose=verbose)
        fig = build_station_map(
            data,
            StationMapOptions(
                backend="matplotlib",
                title=title,
                show_labels=label,
            ),
        )
    except ImportError as exc:
        click.echo(str(exc), err=True)
        sys.exit(1)
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error building station map: {exc}", err=True)
        sys.exit(1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    click.echo(f"Saved station map for {len(rows)} station(s): {output_path}")

    if show:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        plt.show()
    else:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        plt.close(fig)
