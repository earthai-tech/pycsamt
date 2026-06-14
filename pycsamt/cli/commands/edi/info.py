# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt edi info — metadata summary for one EDI file or a whole directory."""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import format_option, no_color_option, verbose_option
from ._base import edi, _get_edi, _get_collection


@edi.command("info")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_OR_DIR",
)
@click.option(
    "--station", "-s",
    default=None,
    metavar="NAME",
    help="Show detail for a single named station (directory mode only).",
)
@click.option(
    "--top", "-n",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Limit output to the first N stations.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def info(
    ctx: click.Context,
    source: Path,
    station: str | None,
    top: int | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Metadata summary for one EDI file or a whole directory.

    In directory mode all EDI files are loaded and a per-station table
    is produced (station, n_freq, freq range, has_tipper, lat, lon, elev).
    In file mode a detailed block-level inspection is shown.

    \b
    Examples:
      pycsamt edi info data/AMT/WILLY_DATA/L18PLT/
      pycsamt edi info data/3edis/ --format csv
      pycsamt edi info data/3edis/ --station HBH03
      pycsamt edi info data/3edis/HBH03.edi
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    if source.is_file():
        _info_single(source, output_format, verbose)
    else:
        _info_collection(source, station, top, output_format, verbose)


# ---------------------------------------------------------------------------
# Single-file display
# ---------------------------------------------------------------------------

def _info_single(path: Path, fmt: str, verbose: int) -> None:
    edi = _get_edi(path, verbose=verbose)
    freq = getattr(edi.Z, "freq", None)
    n_freq = int(getattr(edi.Z, "n_freq", 0) or 0)
    freq_min = float(freq.min()) if freq is not None and freq.size else float("nan")
    freq_max = float(freq.max()) if freq is not None and freq.size else float("nan")

    head = getattr(edi, "Head", None) or {}
    lat  = float(getattr(head, "lat",  None) or float("nan"))
    lon  = float(getattr(head, "lon",  None) or float("nan"))
    elev = float(getattr(head, "elev", None) or float("nan"))

    row = {
        "station":    edi.station or path.stem,
        "path":       str(path),
        "n_freq":     n_freq,
        "freq_min":   freq_min,
        "freq_max":   freq_max,
        "has_tipper": edi.has_tipper,
        "dtype":      edi.dtype,
        "channels":   edi.channels,
        "lat":        lat,
        "lon":        lon,
        "elev":       elev,
    }

    if fmt == "json":
        click.echo(json.dumps(row, indent=2, default=str))
    elif fmt == "csv":
        keys = ["station", "n_freq", "freq_min", "freq_max",
                "has_tipper", "lat", "lon", "elev", "path"]
        click.echo(",".join(keys))
        click.echo(",".join(str(row[k]) for k in keys))
    else:
        click.echo(f"Station     : {row['station']}")
        click.echo(f"File        : {row['path']}")
        click.echo(f"n_freq      : {n_freq}")
        if n_freq > 0:
            click.echo(f"Freq range  : {freq_min:.4g} – {freq_max:.4g} Hz")
        click.echo(f"Has tipper  : {edi.has_tipper}")
        click.echo(f"Data type   : {edi.dtype or '(auto)'}")
        if edi.channels:
            click.echo(f"Channels    : {', '.join(edi.channels)}")
        if not any(c != c for c in [lat, lon]):  # not NaN
            click.echo(f"Lat / Lon   : {lat:.6f}  /  {lon:.6f}")
            click.echo(f"Elevation   : {elev:.1f} m")


# ---------------------------------------------------------------------------
# Collection display
# ---------------------------------------------------------------------------

def _info_collection(
    path: Path,
    station: str | None,
    top: int | None,
    fmt: str,
    verbose: int,
) -> None:
    coll = _get_collection(path, verbose=verbose)

    if len(coll) == 0:
        click.echo(f"No EDI files found in {path}.", err=True)
        return

    rows = []
    for edi in coll:
        freq  = getattr(edi.Z, "freq", None)
        n_freq = int(getattr(edi.Z, "n_freq", 0) or 0)
        freq_min = float(freq.min()) if freq is not None and freq.size else float("nan")
        freq_max = float(freq.max()) if freq is not None and freq.size else float("nan")

        head = getattr(edi, "Head", None)
        lat  = float(getattr(head, "lat",  None) or float("nan")) if head else float("nan")
        lon  = float(getattr(head, "lon",  None) or float("nan")) if head else float("nan")
        elev = float(getattr(head, "elev", None) or float("nan")) if head else float("nan")

        rows.append({
            "station":    edi.station or "?",
            "n_freq":     n_freq,
            "freq_min":   round(freq_min, 4),
            "freq_max":   round(freq_max, 4),
            "has_tipper": edi.has_tipper,
            "lat":        lat,
            "lon":        lon,
            "elev":       elev,
            "path":       edi.path_str,
        })

    if station is not None:
        rows = [r for r in rows if r["station"].upper() == station.upper()]
        if not rows:
            click.echo(
                f"Station {station!r} not found.  "
                f"Available: {', '.join(r['station'] for r in rows[:10])}",
                err=True,
            )
            return

    if top is not None:
        rows = rows[:top]

    if fmt == "json":
        click.echo(json.dumps(rows, indent=2, default=str))
        return

    keys = ["station", "n_freq", "freq_min", "freq_max",
            "has_tipper", "lat", "lon", "elev"]

    if fmt == "csv":
        click.echo(",".join(keys + ["path"]))
        for r in rows:
            click.echo(",".join(str(r[k]) for k in keys + ["path"]))
        return

    # Text table
    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table import Table     # noqa: PLC0415
        tbl = Table(
            "Station", "n_freq", "Freq min", "Freq max",
            "Tipper", "Lat", "Lon", "Elev",
            title=f"EDI collection — {path}  ({len(rows)} stations)",
            show_header=True,
        )
        for r in rows:
            tbl.add_row(
                r["station"],
                str(r["n_freq"]),
                f"{r['freq_min']:.3g}",
                f"{r['freq_max']:.3g}",
                "yes" if r["has_tipper"] else "no",
                f"{r['lat']:.5f}" if r["lat"] == r["lat"] else "—",
                f"{r['lon']:.5f}" if r["lon"] == r["lon"] else "—",
                f"{r['elev']:.0f}" if r["elev"] == r["elev"] else "—",
            )
        Console().print(tbl)
    except ImportError:
        hdr = f"{'Station':<20} {'n_freq':>6} {'freq_min':>10} {'freq_max':>10} " \
              f"{'Tipper':<7} {'Lat':>12} {'Lon':>12} {'Elev':>8}"
        click.echo(hdr)
        click.echo("-" * len(hdr))
        for r in rows:
            click.echo(
                f"{r['station']:<20} {r['n_freq']:>6} "
                f"{r['freq_min']:>10.3g} {r['freq_max']:>10.3g} "
                f"{'yes' if r['has_tipper'] else 'no':<7} "
                f"{r['lat']:>12.5f} {r['lon']:>12.5f} {r['elev']:>8.0f}"
            )
