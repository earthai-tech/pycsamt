# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt avg stations — station positions from AVG and companion .stn file."""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import format_option, no_color_option, verbose_option
from ._base import avg, _get_avg, _load_raw


@avg.command("stations")
@click.argument(
    "source",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    metavar="AVG_FILE",
)
@click.option(
    "--stn-file",
    "stn_file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    metavar="STN_FILE",
    help=(
        "Companion .stn file with UTM / lat-lon coordinates.  "
        "When provided, easting, northing, elevation, and azimuth "
        "are included in the table."
    ),
)
@click.option(
    "--sort-by",
    type=click.Choice(["name", "position", "elev"], case_sensitive=False),
    default="position",
    show_default=True,
    help="Column to sort by.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def stations(
    ctx: click.Context,
    source: Path,
    stn_file: Path | None,
    sort_by: str,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Display the station table for a Zonge AVG file.

    Without ``--stn-file`` only the survey-line positions and station
    names embedded in the AVG file are shown.  With ``--stn-file``
    the full geographic coordinate table (easting, northing, elevation,
    lat/lon) is included.

    \b
    Examples:
      pycsamt avg stations data/avg/K2.AVG
      pycsamt avg stations data/avg/K2.AVG --stn-file data/avg/K2.stn
      pycsamt avg stations data/avg/K2.AVG --stn-file data/avg/K2.stn --format csv
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    try:
        obj = _get_avg(source, verbose=verbose)
    except click.UsageError:
        df, _, _ = _load_raw(source)
        _stations_raw(df, output_format)
        return

    # Basic station info from AVG
    st = obj.info.station
    names  = list(st.names)  if st.names  else []
    values = list(st.values) if st.values is not None else []
    unit   = getattr(st, "unit", "m")

    rows = [
        {"name": n, "position": v, "unit": unit}
        for n, v in zip(names, values)
    ]

    # Attach topography if stn file provided
    topo_attached = False
    topo_frame = None
    if stn_file is not None:
        try:
            obj.add_topography(str(stn_file))
            topo = obj.topo
            if topo is not None and hasattr(topo, "frame") and topo.frame is not None:
                topo_frame = topo.frame
                topo_attached = True
        except Exception as exc:  # noqa: BLE001
            click.echo(f"Warning: could not load {stn_file.name}: {exc}", err=True)

    if topo_attached and topo_frame is not None:
        import numpy as np  # noqa: PLC0415
        tf = topo_frame.copy()
        stn_col = (
            tf["station"].to_numpy(float, na_value=float("nan"))
            if "station" in tf.columns
            else np.arange(len(tf), dtype=float)
        )
        topo_cols = [c for c in ("easting", "northing", "elevation",
                                  "latitude", "longitude") if c in tf.columns]
        for r in rows:
            pos = r.get("position")
            if pos is not None:
                # nearest-neighbour match by position value
                idx = int(np.argmin(np.abs(stn_col - float(pos))))
                row_topo = tf.iloc[idx]
            else:
                # fall back to first row
                row_topo = tf.iloc[0]
            for col in topo_cols:
                r[col] = float(row_topo[col])

    # Sort
    sort_map = {"name": "name", "position": "position", "elev": "elevation"}
    key = sort_map.get(sort_by, "position")
    rows.sort(key=lambda r: r.get(key, 0) or 0)

    if output_format == "json":
        click.echo(json.dumps(rows, indent=2, default=str))
        return

    has_topo = topo_attached and any("easting" in r for r in rows)
    base_cols = ["name", "position"]
    topo_cols = ["easting", "northing", "elevation", "latitude", "longitude"]
    display_cols = base_cols + (topo_cols if has_topo else [])

    if output_format == "csv":
        click.echo(",".join(display_cols))
        for r in rows:
            click.echo(",".join(str(r.get(c, "")) for c in display_cols))
        return

    # Text table
    click.echo(f"\nStations — {source.name}  (unit: {unit})\n")
    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table import Table     # noqa: PLC0415
        tbl = Table("Name", f"Position ({unit})")
        if has_topo:
            for c in ["Easting", "Northing", "Elev (m)", "Lat", "Lon"]:
                tbl.add_column(c, justify="right")
        for r in rows:
            row_vals = [str(r.get("name", "?")), f"{r.get('position', 0):.1f}"]
            if has_topo:
                row_vals += [
                    f"{r.get('easting',  float('nan')):.1f}",
                    f"{r.get('northing', float('nan')):.1f}",
                    f"{r.get('elevation',float('nan')):.1f}",
                    f"{r.get('latitude', float('nan')):.5f}",
                    f"{r.get('longitude',float('nan')):.5f}",
                ]
            tbl.add_row(*row_vals)
        Console().print(tbl)
    except ImportError:
        hdr = f"{'Name':<10} {'Pos(m)':>10}"
        if has_topo:
            hdr += f" {'E(m)':>12} {'N(m)':>12} {'Elev(m)':>9} {'Lat':>12} {'Lon':>12}"
        click.echo(hdr)
        click.echo("-" * len(hdr))
        for r in rows:
            line = f"{str(r.get('name','?')):<10} {r.get('position', 0):>10.1f}"
            if has_topo:
                line += (
                    f" {r.get('easting',  float('nan')):>12.1f}"
                    f" {r.get('northing', float('nan')):>12.1f}"
                    f" {r.get('elevation',float('nan')):>9.1f}"
                    f" {r.get('latitude', float('nan')):>12.5f}"
                    f" {r.get('longitude',float('nan')):>12.5f}"
                )
            click.echo(line)

    click.echo(f"\n{len(rows)} station(s)")


def _stations_raw(df, fmt):
    """Minimal station table from raw DataFrame."""
    if "station" not in df.columns:
        click.echo("No station column in data.", err=True)
        return
    sta = sorted(df["station"].unique())
    rows = [{"station": str(s)} for s in sta]
    if fmt == "json":
        click.echo(json.dumps(rows, indent=2))
    elif fmt == "csv":
        click.echo("station")
        for r in rows:
            click.echo(r["station"])
    else:
        for r in rows:
            click.echo(r["station"])
