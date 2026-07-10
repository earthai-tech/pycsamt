# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt jones stations — station coordinate and metadata table."""

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
from ._base import _get_collection, jones


@jones.command("stations")
@click.argument(
    "source",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="J_DIR",
)
@click.option(
    "--sort-by",
    type=click.Choice(["station", "lat", "lon", "elev", "az", "n_freq"],
                      case_sensitive=False),
    default="station",
    show_default=True,
    help="Column to sort the output table by.",
)
@click.option(
    "--top", "-n",
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
    sort_by: str,
    top: int | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Display station metadata for all J-files in a directory.

    Columns: station, lat, lon, elev (m), azimuth (°), n_freq,
    has_z (impedance), has_r (rho/phase), has_t (tipper).

    \b
    Examples:
      pycsamt jones stations data/j/
      pycsamt jones stations data/j/ --sort-by lat --format csv
      pycsamt jones stations data/j/ --sort-by n_freq --top 5
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    coll = _get_collection(source, verbose=verbose)
    if len(coll) == 0:
        click.echo(f"No J-files found in {source}.", err=True)
        return

    rows = coll.summary()

    # Sort
    _key = {
        "station": lambda r: str(r.get("station", "")),
        "lat":     lambda r: float(r.get("lat",    0) or 0),
        "lon":     lambda r: float(r.get("lon",    0) or 0),
        "elev":    lambda r: float(r.get("elev",   0) or 0),
        "az":      lambda r: float(r.get("az",     0) or 0),
        "n_freq":  lambda r: int(  r.get("n_freq", 0) or 0),
    }
    rows.sort(key=_key[sort_by])

    if top is not None:
        rows = rows[:top]

    if output_format == "json":
        click.echo(json.dumps(rows, indent=2, default=str))
        return

    display_cols = [
        "station", "lat", "lon", "az", "n_freq", "has_z", "has_r", "has_t"
    ]

    if output_format == "csv":
        click.echo(",".join(display_cols))
        for r in rows:
            vals = []
            for c in display_cols:
                v = r.get(c, "")
                vals.append(str(v))
            click.echo(",".join(vals))
        return

    # Rich text
    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table import Table  # noqa: PLC0415
        tbl = Table(
            title=f"Station metadata — {source}  ({len(rows)} stations)",
        )
        labels = {
            "station": "Station", "lat": "Lat",  "lon": "Lon",
            "az": "Az(°)", "n_freq": "n_freq",
            "has_z": "Z", "has_r": "R/φ", "has_t": "T",
        }
        for c in display_cols:
            tbl.add_column(labels.get(c, c))
        for r in rows:
            row_vals = []
            for c in display_cols:
                v = r.get(c, "")
                if isinstance(v, float) and c in ("lat", "lon"):
                    row_vals.append(f"{v:.5f}")
                elif isinstance(v, bool):
                    row_vals.append("yes" if v else "no")
                else:
                    row_vals.append(str(v) if v is not None else "—")
            tbl.add_row(*row_vals)
        Console().print(tbl)
    except ImportError:
        hdr = (
            f"{'Station':<18} {'Lat':>12} {'Lon':>12} {'Az':>8} "
            f"{'n_freq':>6} {'Z':<4} {'R/φ':<4} {'T':<4}"
        )
        click.echo(hdr)
        click.echo("-" * len(hdr))
        for r in rows:
            lat = r.get("lat", float("nan")) or float("nan")
            lon = r.get("lon", float("nan")) or float("nan")
            az  = r.get("az",  float("nan")) or float("nan")
            click.echo(
                f"{str(r.get('station', '?')):<18} "
                f"{lat:>12.5f} {lon:>12.5f} {az:>8.1f} "
                f"{r.get('n_freq', 0):>6} "
                f"{'y' if r.get('has_z') else 'n':<4}"
                f"{'y' if r.get('has_r') else 'n':<4}"
                f"{'y' if r.get('has_t') else 'n':<4}"
            )
    click.echo(f"\n{len(rows)} station(s) total")
