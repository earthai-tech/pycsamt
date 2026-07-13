# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt edi stations — station coordinate table."""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    no_color_option,
    verbose_option,
)
from ._base import _get_collection, edi


@edi.command("stations")
@click.argument(
    "source",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="EDI_DIR",
)
@click.option(
    "--pattern",
    "-p",
    default=None,
    metavar="GLOB",
    help="Filter station names by a glob pattern (e.g. 'S0*', '18-02*').",
)
@click.option(
    "--sort-by",
    type=click.Choice(
        ["station", "lat", "lon", "elev"], case_sensitive=False
    ),
    default="station",
    show_default=True,
    help="Column to sort the output table by.",
)
@click.option(
    "--top",
    "-n",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Limit output to the first N rows.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def stations(
    ctx: click.Context,
    source: Path,
    pattern: str | None,
    sort_by: str,
    top: int | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Display a station coordinate table for all EDIs in a directory.

    Columns: station, lat, lon, elev (m), easting (m), northing (m),
    UTM zone, and file path.

    \b
    Examples:
      pycsamt edi stations data/AMT/WILLY_DATA/L18PLT/
      pycsamt edi stations data/AMT/ --sort-by lat --format csv
      pycsamt edi stations data/AMT/ --pattern "18-02*" --top 10
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    import fnmatch  # noqa: PLC0415

    from pycsamt.seg.survey import Stations  # noqa: PLC0415

    coll = _get_collection(source, verbose=verbose)
    if len(coll) == 0:
        click.echo(f"No EDI files found in {source}.", err=True)
        return

    st = Stations(coll)
    rows = [r for r in st.table() if "ed" not in r or True]

    # strip the heavy 'ed' key (EDIFile object) for display
    rows = [{k: v for k, v in r.items() if k != "ed"} for r in rows]

    # Filter by pattern
    if pattern is not None:
        rows = [
            r
            for r in rows
            if fnmatch.fnmatch(str(r.get("station", "")), pattern)
        ]

    # Sort
    _sort_key = {
        "station": lambda r: str(r.get("station", "")),
        "lat": lambda r: float(r.get("lat", 0) or 0),
        "lon": lambda r: float(r.get("lon", 0) or 0),
        "elev": lambda r: float(r.get("elev", 0) or 0),
    }
    rows.sort(key=_sort_key[sort_by])

    if top is not None:
        rows = rows[:top]

    if output_format == "json":
        click.echo(json.dumps(rows, indent=2, default=str))
        return

    display_cols = ["station", "lat", "lon", "elev", "e", "n", "zone", "path"]
    display_cols = [c for c in display_cols if any(c in r for r in rows)]

    if output_format == "csv":
        click.echo(",".join(display_cols))
        for r in rows:
            click.echo(",".join(str(r.get(c, "")) for c in display_cols))
        return

    # Rich text table
    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table import Table  # noqa: PLC0415

        tbl = Table(
            title=f"Station coordinates — {source}  ({len(rows)} stations)",
            show_header=True,
        )
        text_cols = ["station", "lat", "lon", "elev", "zone"]
        for c in text_cols:
            if any(c in r for r in rows):
                tbl.add_column(c)
        for r in rows:
            tbl.add_row(
                *[
                    f"{r.get(c, ''):.5f}"
                    if isinstance(r.get(c), float)
                    else str(r.get(c, ""))
                    for c in text_cols
                    if any(c in row for row in rows)
                ]
            )
        Console().print(tbl)
    except ImportError:
        hdr = (
            f"{'Station':<22} {'Lat':>12} {'Lon':>12} {'Elev':>8} {'Zone':<6}"
        )
        click.echo(hdr)
        click.echo("-" * len(hdr))
        for r in rows:
            lat = r.get("lat", float("nan"))
            lon = r.get("lon", float("nan"))
            elev = r.get("elev", float("nan"))
            zone = r.get("zone", "")
            click.echo(
                f"{str(r.get('station', '?')):<22} "
                f"{lat:>12.5f} {lon:>12.5f} {elev:>8.1f} {zone:<6}"
            )

    click.echo(f"\n{len(rows)} station(s) total")
