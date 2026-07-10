# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt edi profile — survey profile geometry (bearing, step, distances)."""

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


@edi.command("profile")
@click.argument(
    "source",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="EDI_DIR",
)
@click.option(
    "--bearing-method",
    type=click.Choice(["endpoints", "linear"], case_sensitive=False),
    default="endpoints",
    show_default=True,
    help=(
        "Method for estimating the profile bearing.  "
        "'endpoints' uses the first-to-last station vector.  "
        "'linear' fits the best axis by SVD of projected coordinates."
    ),
)
@click.option(
    "--step-method",
    type=click.Choice(["mean", "median"], case_sensitive=False),
    default="mean",
    show_default=True,
    help="Aggregation for the inter-station spacing scalar.",
)
@click.option(
    "--distances",
    is_flag=True,
    default=False,
    help="Include per-station cumulative along-profile distances in output.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def profile(
    ctx: click.Context,
    source: Path,
    bearing_method: str,
    step_method: str,
    distances: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Survey profile geometry: bearing, spacing, and station offsets.

    Builds an :class:`~pycsamt.seg.survey.EDIProfile` from the collection
    and reports the survey azimuth, mean inter-station step, and optionally
    the per-station cumulative distance along the profile.

    \b
    Examples:
      pycsamt edi profile data/AMT/WILLY_DATA/L18PLT/
      pycsamt edi profile data/AMT/WILLY_DATA/L18PLT/ --bearing-method linear
      pycsamt edi profile data/AMT/WILLY_DATA/L18PLT/ --distances --format csv
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from pycsamt.seg.survey import EDIProfile  # noqa: PLC0415

    coll = _get_collection(source, verbose=verbose)
    if len(coll) == 0:
        click.echo(f"No EDI files found in {source}.", err=True)
        return

    prof = EDIProfile(coll)
    bearing  = prof.get_bearing(method=bearing_method)
    step_m   = prof.get_step(method=step_method)
    n_sta    = len(prof.stations)
    tbl      = prof.as_table()

    summary = {
        "source":          str(source),
        "n_stations":      n_sta,
        "bearing_deg":     round(bearing, 2) if bearing is not None else None,
        "bearing_method":  bearing_method,
        "step_m":          round(float(step_m), 1) if step_m is not None else None,
        "step_method":     step_method,
    }

    if distances:
        dist = prof.distance if hasattr(prof, "distance") else None
        dist_list = dist.tolist() if dist is not None else None
        summary["distances_m"] = dist_list
        summary["station_table"] = [
            {k: v for k, v in r.items() if k not in ("easting", "northing")}
            for r in tbl
        ]

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2, default=str))
        return

    if output_format == "csv":
        rows = tbl if tbl else []
        if rows:
            cols = [c for c in rows[0].keys() if c not in ("easting", "northing")]
            click.echo(",".join(cols))
            for r in rows:
                click.echo(",".join(str(r.get(c, "")) for c in cols))
        return

    # Text / rich output
    click.echo(f"\nSurvey profile — {source}")
    click.echo(f"  Stations   : {n_sta}")
    if bearing is not None:
        click.echo(f"  Bearing    : {bearing:.1f}°  ({bearing_method})")
    else:
        click.echo("  Bearing    : (insufficient data)")
    if step_m is not None:
        click.echo(f"  Step (avg) : {float(step_m):.0f} m  ({step_method})")

    if distances and tbl:
        dist = prof.distance if hasattr(prof, "distance") else None
        click.echo()
        hdr = f"  {'Station':<22} {'Lat':>10} {'Lon':>10} {'Elev':>7} {'Dist (m)':>10}"
        click.echo(hdr)
        click.echo("  " + "-" * (len(hdr) - 2))
        for i, r in enumerate(tbl):
            d = f"{dist[i]:>10.0f}" if dist is not None and i < len(dist) else f"{'—':>10}"
            click.echo(
                f"  {str(r.get('station', '?')):<22} "
                f"{r.get('lat', float('nan')):>10.5f} "
                f"{r.get('lon', float('nan')):>10.5f} "
                f"{r.get('elev', float('nan')):>7.1f} "
                f"{d}"
            )
