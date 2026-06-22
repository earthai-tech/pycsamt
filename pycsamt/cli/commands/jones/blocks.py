# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt jones blocks — inspect data blocks inside a J-file."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import format_option, no_color_option, verbose_option
from ._base import jones, _get_jfile


@jones.command("blocks")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="J_FILE",
)
@click.option(
    "--kind", "-k",
    default=None,
    metavar="KIND",
    help=(
        "Filter blocks by data-type kind.  "
        "R = rho/phase,  Z = impedance,  T = tipper/GDS,  "
        "C = Schmucker C.  Case-insensitive."
    ),
)
@click.option(
    "--comp", "-c",
    default=None,
    metavar="COMP",
    help=(
        "Filter blocks by tensor component, e.g. XY, YX, ZX, ZY.  "
        "Case-insensitive."
    ),
)
@click.option(
    "--qa",
    is_flag=True,
    default=False,
    help="Show data-quality summary for each block (rejection rate, period range).",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def blocks(
    ctx: click.Context,
    source: Path,
    kind: str | None,
    comp: str | None,
    qa: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Inspect the data blocks contained in a single J-file.

    Reports each block's kind (R/Z/T/C), tensor component (XX/XY/YX/YY
    for impedance, ZX/ZY for tipper), number of rows, and period range.
    Use ``--qa`` to add data-quality indicators.

    \b
    Examples:
      pycsamt jones blocks data/j/S01.j
      pycsamt jones blocks data/j/S01.j --kind R
      pycsamt jones blocks data/j/S01.j --comp XY --qa
      pycsamt jones blocks data/j/S01.j --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    jf = _get_jfile(source, verbose=verbose)
    blks = jf.blocks

    if blks is None or blks.n == 0:
        click.echo(f"No data blocks found in {source}.", err=True)
        return

    # Filter
    selected = blks.blocks
    if kind is not None:
        selected = [b for b in selected
                    if (getattr(b, "kind", "") or "").upper() == kind.upper()]
    if comp is not None:
        selected = [b for b in selected
                    if (getattr(b, "comp", "") or "").upper() == comp.upper()]

    if not selected:
        click.echo(
            f"No blocks match kind={kind!r} comp={comp!r}.", err=True
        )
        sys.exit(1)

    rows: list[dict] = []
    for b in selected:
        bkind   = getattr(b, "kind",  "?") or "?"
        bcomp   = getattr(b, "comp",  "?") or "?"
        nrows   = int(getattr(b, "nrows", 0) or 0)
        periods = getattr(b, "periods", None)
        p_min   = float(periods.min()) if periods is not None and periods.size else float("nan")
        p_max   = float(periods.max()) if periods is not None and periods.size else float("nan")

        row: dict = {
            "kind":       bkind,
            "comp":       bcomp,
            "token":      f"{bkind}{bcomp}",
            "nrows":      nrows,
            "period_min": round(p_min, 6),
            "period_max": round(p_max, 6),
        }

        if qa:
            try:
                qa_info = b.qa_summary()
                row["qa"] = qa_info
            except Exception:  # noqa: BLE001
                row["qa"] = {}

        rows.append(row)

    station = jf.station or source.stem

    if output_format == "json":
        click.echo(json.dumps({
            "station": station,
            "source":  str(source),
            "n_blocks": len(rows),
            "blocks":  rows,
        }, indent=2, default=str))
        return

    if output_format == "csv":
        csv_cols = ["token", "kind", "comp", "nrows", "period_min", "period_max"]
        click.echo(",".join(csv_cols))
        for r in rows:
            click.echo(",".join(str(r.get(c, "")) for c in csv_cols))
        return

    # Text
    click.echo(f"\nStation : {station}  ({source})")
    click.echo(f"Blocks  : {len(rows)}")
    click.echo()
    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table import Table     # noqa: PLC0415
        tbl = Table("Token", "Kind", "Comp", "Rows", "T_min (s)", "T_max (s)")
        if qa:
            tbl.add_column("QA")
        for r in rows:
            row_vals = [
                r["token"], r["kind"], r["comp"], str(r["nrows"]),
                f"{r['period_min']:.4g}", f"{r['period_max']:.4g}",
            ]
            if qa:
                qa_info = r.get("qa", {})
                rej = qa_info.get("rejected", 0)
                n   = qa_info.get("n", r["nrows"])
                row_vals.append(f"{rej}/{n} rej")
            tbl.add_row(*row_vals)
        Console().print(tbl)
    except ImportError:
        hdr = (
            f"{'Token':<8} {'Kind':<6} {'Comp':<6} "
            f"{'Rows':>5} {'T_min(s)':>12} {'T_max(s)':>12}"
        )
        if qa:
            hdr += f"  {'QA'}"
        click.echo(hdr)
        click.echo("-" * len(hdr))
        for r in rows:
            line = (
                f"{r['token']:<8} {r['kind']:<6} {r['comp']:<6} "
                f"{r['nrows']:>5} {r['period_min']:>12.4g} {r['period_max']:>12.4g}"
            )
            if qa:
                qa_info = r.get("qa", {})
                rej = qa_info.get("rejected", 0)
                n   = qa_info.get("n", r["nrows"])
                line += f"  {rej}/{n} rej"
            click.echo(line)
