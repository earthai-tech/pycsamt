# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt avg info — metadata summary for a Zonge AVG file."""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import format_option, no_color_option, verbose_option
from ._base import avg, _get_avg, _load_raw


@avg.command("info")
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
    help="Optional companion .stn file to attach station coordinates.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def info(
    ctx: click.Context,
    source: Path,
    stn_file: Path | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Display metadata for a Zonge AVG file.

    Reports project name, survey type, hardware, data format,
    station coverage, frequency range, and component availability.

    \b
    Examples:
      pycsamt avg info data/avg/K2.AVG
      pycsamt avg info data/avg/K2.AVG --stn-file data/avg/K2.stn
      pycsamt avg info data/avg/K2.AVG --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    try:
        obj = _get_avg(source, verbose=verbose)
        _info_rich(obj, stn_file, source, output_format)
    except click.UsageError:
        # Legacy file without xarray — fall back to raw loader
        df, meta, kind = _load_raw(source)
        _info_raw(df, meta, kind, source, output_format)


def _info_rich(obj, stn_file, source, fmt):
    """Display info using the full AVG object (K2 path)."""
    if stn_file is not None:
        try:
            obj.add_topography(str(stn_file))
        except Exception:  # noqa: BLE001
            pass

    s = obj.summary
    hdr = obj.info.header if obj.info else None
    ann = getattr(hdr, "annotation", None) if hdr else None
    cfg = getattr(hdr, "config",     None) if hdr else None
    hw  = getattr(hdr, "hardware",   None) if hdr else None

    comps = sorted(obj.df["comp"].unique()) if obj.df is not None else []

    row = {
        "source":         str(source),
        "data_kind":      s.data_kind,
        "project":        s.project,
        "survey_type":    getattr(cfg, "survey_type", None) or s.survey_type,
        "line_name":      getattr(ann, "line_name",   None) or getattr(cfg, "line_name", None),
        "operator":       getattr(ann, "operator",    None),
        "n_stations":     s.num_stations,
        "n_frequencies":  s.num_frequencies,
        "station_range":  s.station_range,
        "frequency_range": s.frequency_range,
        "components":     comps,
        "total_rows":     s.total_rows,
        "hardware":       str(hw) if hw else None,
    }

    if fmt == "json":
        click.echo(json.dumps(row, indent=2, default=str))
        return

    if fmt == "csv":
        keys = ["project", "survey_type", "n_stations", "n_frequencies",
                "station_range", "frequency_range", "data_kind"]
        click.echo(",".join(keys))
        click.echo(",".join(str(row.get(k, "")) for k in keys))
        return

    click.echo(f"File          : {source}")
    click.echo(f"Format        : {row['data_kind']}")
    click.echo(f"Project       : {row['project']}")
    if row["survey_type"]:
        click.echo(f"Survey type   : {row['survey_type']}")
    if row["line_name"]:
        click.echo(f"Line          : {row['line_name']}")
    if row["operator"]:
        click.echo(f"Operator      : {row['operator']}")
    click.echo(f"Stations      : {row['n_stations']}  ({row['station_range']})")
    click.echo(f"Frequencies   : {row['n_frequencies']}  ({row['frequency_range']})")
    if comps:
        click.echo(f"Components    : {', '.join(comps)}")
    click.echo(f"Total rows    : {row['total_rows']}")


def _info_raw(df, meta, kind, source, fmt):
    """Display info from raw DataFrame (K1 legacy path without xarray)."""
    n_sta  = int(df["station"].nunique()) if "station" in df.columns else 0
    n_freq = int(df["freq"].nunique())    if "freq"    in df.columns else 0
    comps  = sorted(df["comp"].unique())  if "comp"    in df.columns else []
    freq   = df["freq"] if "freq" in df.columns else None
    f_min  = float(freq.min()) if freq is not None else float("nan")
    f_max  = float(freq.max()) if freq is not None else float("nan")

    row = {
        "source":        str(source),
        "data_kind":     f"Kind-{kind} (legacy)",
        "n_stations":    n_sta,
        "n_frequencies": n_freq,
        "frequency_range": f"{f_min:.4g} – {f_max:.4g} Hz" if freq is not None else "N/A",
        "components":    comps,
        "total_rows":    len(df),
        "note":          "xarray not installed — limited metadata available for legacy files",
    }

    if fmt == "json":
        click.echo(json.dumps(row, indent=2, default=str))
        return

    click.echo(f"File          : {source}")
    click.echo(f"Format        : {row['data_kind']}")
    click.echo(f"Stations      : {n_sta}")
    click.echo(f"Frequencies   : {n_freq}  ({row['frequency_range']})")
    if comps:
        click.echo(f"Components    : {', '.join(comps)}")
    click.echo(f"Total rows    : {len(df)}")
    click.echo(f"Note          : {row['note']}")
