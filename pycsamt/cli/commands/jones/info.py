# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt jones info — metadata summary for one J-file or a directory."""

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
from ._base import _get_collection, _get_jfile, jones


@jones.command("info")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="J_OR_DIR",
)
@click.option(
    "--top",
    "-n",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Limit table to the first N stations (directory mode).",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def info(
    ctx: click.Context,
    source: Path,
    top: int | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Metadata summary for one Jones J-file or a whole directory.

    Single-file mode shows station, n_freq, frequency range, data
    blocks, coordinates, and banner provenance.
    Directory mode produces a per-station table.

    \b
    Examples:
      pycsamt jones info data/j/S01.j
      pycsamt jones info data/j/
      pycsamt jones info data/j/ --format csv
      pycsamt jones info data/j/ --top 10
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    if source.is_file():
        _info_single(source, output_format, verbose)
    else:
        _info_collection(source, top, output_format, verbose)


# ---------------------------------------------------------------------------
# Single-file display
# ---------------------------------------------------------------------------


def _info_single(path: Path, fmt: str, verbose: int) -> None:
    jf = _get_jfile(path, verbose=verbose)

    freq = jf.freq
    n_freq = int(jf.n_freq or 0)
    freq_min = (
        float(freq.min()) if freq is not None and freq.size else float("nan")
    )
    freq_max = (
        float(freq.max()) if freq is not None and freq.size else float("nan")
    )

    banner = getattr(getattr(jf, "heads", None), "banner", None)
    software = getattr(banner, "software", None) if banner else None
    bdate = getattr(banner, "date", None) if banner else None

    blocks = jf.blocks
    block_kinds = []
    if blocks is not None:
        for b in blocks.blocks:
            kind = getattr(b, "kind", None) or "?"
            comp = getattr(b, "comp", None) or "?"
            nr = int(getattr(b, "nrows", 0) or 0)
            block_kinds.append({"kind": kind, "comp": comp, "nrows": nr})

    row = {
        "station": jf.station or path.stem,
        "path": str(path),
        "n_freq": n_freq,
        "freq_min": round(freq_min, 4),
        "freq_max": round(freq_max, 4),
        "has_z": jf.Z is not None and getattr(jf.Z, "z", None) is not None,
        "has_r": jf.Res is not None,
        "has_t": jf.Tip is not None,
        "lat": jf.lat,
        "lon": jf.lon,
        "elev": jf.elev,
        "azimuth": jf.azimuth,
        "software": software,
        "date": bdate,
        "blocks": block_kinds,
    }

    if fmt == "json":
        click.echo(json.dumps(row, indent=2, default=str))
        return

    if fmt == "csv":
        keys = [
            "station",
            "n_freq",
            "freq_min",
            "freq_max",
            "has_z",
            "has_r",
            "has_t",
            "lat",
            "lon",
            "elev",
            "path",
        ]
        click.echo(",".join(keys))
        click.echo(",".join(str(row.get(k, "")) for k in keys))
        return

    click.echo(f"Station   : {row['station']}")
    click.echo(f"File      : {row['path']}")
    click.echo(f"n_freq    : {n_freq}")
    if n_freq > 0:
        click.echo(f"Freq range: {freq_min:.4g} – {freq_max:.4g} Hz")
    click.echo(
        f"Data      : Z={row['has_z']}  R/φ={row['has_r']}  T={row['has_t']}"
    )
    if row["lat"] is not None:
        click.echo(f"Lat / Lon : {row['lat']:.6f}  /  {row['lon']:.6f}")
    if row["elev"] is not None:
        click.echo(f"Elevation : {row['elev']:.1f} m")
    if row["azimuth"] is not None:
        click.echo(f"Azimuth   : {row['azimuth']:.1f}°")
    if software:
        click.echo(f"Software  : {software}  ({bdate or '?'})")
    if block_kinds:
        click.echo(f"Blocks    : {len(block_kinds)}")
        if verbose >= 1:
            for bk in block_kinds:
                click.echo(f"  {bk['kind']}{bk['comp']}  ({bk['nrows']} rows)")


# ---------------------------------------------------------------------------
# Collection display
# ---------------------------------------------------------------------------


def _info_collection(
    path: Path, top: int | None, fmt: str, verbose: int
) -> None:
    coll = _get_collection(path, verbose=verbose)
    if len(coll) == 0:
        click.echo(f"No J-files found in {path}.", err=True)
        return

    rows = coll.summary()
    if top is not None:
        rows = rows[:top]

    if fmt == "json":
        click.echo(json.dumps(rows, indent=2, default=str))
        return

    keys = [
        "station",
        "n_freq",
        "has_z",
        "has_r",
        "has_t",
        "lat",
        "lon",
        "az",
    ]
    keys = [k for k in keys if any(k in r for r in rows)]

    if fmt == "csv":
        click.echo(",".join(keys))
        for r in rows:
            click.echo(",".join(str(r.get(k, "")) for k in keys))
        return

    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table import Table  # noqa: PLC0415

        tbl = Table(
            title=f"Jones J collection — {path}  ({len(rows)} stations)",
        )
        col_labels = {
            "station": "Station",
            "n_freq": "n_freq",
            "has_z": "Z",
            "has_r": "R/φ",
            "has_t": "T",
            "lat": "Lat",
            "lon": "Lon",
            "az": "Az(°)",
        }
        for k in keys:
            tbl.add_column(col_labels.get(k, k))
        for r in rows:
            tbl.add_row(
                *[
                    f"{r.get(k, ''):.5f}"
                    if isinstance(r.get(k), float) and k in ("lat", "lon")
                    else (
                        "yes"
                        if r.get(k) is True
                        else "no"
                        if r.get(k) is False
                        else str(r.get(k, ""))
                    )
                    for k in keys
                ]
            )
        Console().print(tbl)
    except ImportError:
        hdr = (
            f"{'Station':<18} {'n_freq':>6} {'Z':<4} {'R/φ':<4} {'T':<4} "
            f"{'Lat':>12} {'Lon':>12} {'Az':>8}"
        )
        click.echo(hdr)
        click.echo("-" * len(hdr))
        for r in rows:
            lat = r.get("lat", float("nan"))
            lon = r.get("lon", float("nan"))
            az = r.get("az", float("nan"))
            click.echo(
                f"{str(r.get('station', '?')):<18} "
                f"{r.get('n_freq', 0):>6} "
                f"{'y' if r.get('has_z') else 'n':<4} "
                f"{'y' if r.get('has_r') else 'n':<4} "
                f"{'y' if r.get('has_t') else 'n':<4} "
                f"{lat:>12.5f} {lon:>12.5f} {az:>8.1f}"
            )
