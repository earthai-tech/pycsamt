# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt avg validate — data quality assessment for Zonge AVG files."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import format_option, no_color_option, verbose_option
from ._base import avg, _get_avg, _load_raw

# QC column registry: (col_name, label, unit, threshold_warn)
_QC_COLS = [
    ("pc_emag",  "E-mag % err",   "%",   5.0),
    ("pc_hmag",  "H-mag % err",   "%",   5.0),
    ("pc_rho",   "ρa % err",      "%",  10.0),
    ("s_ephz",   "E-phase σ",     "mrad", 10.0),
    ("s_hphz",   "H-phase σ",     "mrad", 10.0),
    ("s_phz",    "phase σ",       "mrad", 10.0),
    ("z.%err",   "Z % err",       "%",   10.0),
]


@avg.command("validate")
@click.argument(
    "source",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    metavar="AVG_FILE",
)
@click.option(
    "--comp",
    default=None,
    metavar="COMP",
    help="Filter to one component (e.g. ZXY, ZYX, RXY).  Default: all.",
)
@click.option(
    "--threshold",
    "warn_threshold",
    type=float,
    default=None,
    metavar="PCT",
    help=(
        "Override the default QC warning threshold (%).  "
        "Stations whose mean exceeds this are flagged."
    ),
)
@click.option(
    "--top", "-n",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Show only the N worst stations (sorted by mean % error).",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def validate(
    ctx: click.Context,
    source: Path,
    comp: str | None,
    warn_threshold: float | None,
    top: int | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Assess data quality in a Zonge AVG file.

    Computes per-station means for QC metrics already embedded in
    the AVG file: percent error for E-mag, H-mag, and apparent
    resistivity; phase standard deviations for E-field and H-field;
    and the composite Z % error.  Stations exceeding the warning
    threshold are flagged.

    Works on both modern (kind-2) and legacy (kind-1) AVG files.

    \b
    Examples:
      pycsamt avg validate data/avg/K2.AVG
      pycsamt avg validate data/avg/K2.AVG --comp ZXY
      pycsamt avg validate data/avg/K2.AVG --threshold 8.0 --top 10
      pycsamt avg validate data/avg/K2.AVG --format csv
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    # Load the raw DataFrame (works for both K1 and K2)
    try:
        obj = _get_avg(source, verbose=verbose)
        df  = obj.df
    except click.UsageError:
        df, _, _ = _load_raw(source)

    if df is None or df.empty:
        click.echo("No data loaded.", err=True)
        sys.exit(1)

    if comp is not None:
        mask = df["comp"].str.upper() == comp.upper() if "comp" in df.columns else True
        df   = df[mask]
        if df.empty:
            click.echo(f"No rows match component {comp!r}.", err=True)
            sys.exit(1)

    # Find which QC columns are present
    present = [(col, lbl, unit, thr) for col, lbl, unit, thr in _QC_COLS
               if col in df.columns]

    if not present:
        click.echo(
            "No QC columns found in this file (pc_emag, s_ephz, etc. missing).",
            err=True,
        )
        sys.exit(1)

    # Per-station aggregation
    import pandas as pd  # noqa: PLC0415
    import numpy as np   # noqa: PLC0415

    gb  = df.groupby("station") if "station" in df.columns else df.groupby(df.index)
    agg = gb[[c for c, *_ in present]].agg(["mean", "max"]).round(3)
    agg.columns = ["_".join(c) for c in agg.columns]

    rows = []
    for stn, row in agg.iterrows():
        entry = {"station": str(stn)}
        flagged = False
        for col, lbl, unit, default_thr in present:
            thr   = warn_threshold if warn_threshold is not None else default_thr
            mean_ = float(row.get(f"{col}_mean", float("nan")))
            max_  = float(row.get(f"{col}_max",  float("nan")))
            entry[f"{col}_mean"] = mean_
            entry[f"{col}_max"]  = max_
            if mean_ > thr:
                flagged = True
        entry["flagged"] = flagged
        rows.append(entry)

    # Sort worst first, then alpha
    rows.sort(key=lambda r: (
        not r["flagged"],
        -max(
            r.get(f"{col}_mean", 0) or 0
            for col, *_ in present
        ),
    ))

    if top is not None:
        rows = rows[:top]

    n_flagged = sum(1 for r in rows if r["flagged"])

    if output_format == "json":
        click.echo(json.dumps({
            "source":     str(source),
            "n_stations": len(rows),
            "n_flagged":  n_flagged,
            "qc_columns": [col for col, *_ in present],
            "stations":   rows,
        }, indent=2, default=str))
        return

    if output_format == "csv":
        cols = ["station"] + [
            f"{col}_mean" for col, *_ in present
        ] + ["flagged"]
        click.echo(",".join(cols))
        for r in rows:
            click.echo(",".join(str(r.get(c, "")) for c in cols))
        return

    # Rich text table
    click.echo(
        f"\nQC summary — {source.name}  "
        f"({len(rows)} stations, {n_flagged} flagged)\n"
    )
    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table import Table     # noqa: PLC0415
        tbl = Table(show_header=True)
        tbl.add_column("Station", style="cyan")
        for _, lbl, unit, _ in present:
            tbl.add_column(f"{lbl} ({unit})", justify="right")
        tbl.add_column("Flag")
        for r in rows:
            vals = [str(r["station"])]
            for col, lbl, unit, thr in present:
                v = r.get(f"{col}_mean", float("nan"))
                t = warn_threshold if warn_threshold is not None else thr
                color = "red" if v > t else "green"
                vals.append(f"[{color}]{v:.3f}[/{color}]")
            vals.append("[red]✗[/red]" if r["flagged"] else "[green]✓[/green]")
            tbl.add_row(*vals)
        Console().print(tbl)
    except ImportError:
        hdr_cols = ["Station"] + [lbl for _, lbl, *_ in present] + ["OK?"]
        widths   = [12] + [14] * len(present) + [5]
        click.echo("  ".join(f"{h:<{w}}" for h, w in zip(hdr_cols, widths)))
        click.echo("-" * (sum(widths) + 2 * len(widths)))
        for r in rows:
            vals = [str(r["station"])[:12]]
            for col, *_ in present:
                v = r.get(f"{col}_mean", float("nan"))
                vals.append(f"{v:>14.3f}")
            vals.append("  ✗" if r["flagged"] else "  ✓")
            click.echo("  ".join(f"{v:<{w}}" for v, w in zip(vals, widths)))
