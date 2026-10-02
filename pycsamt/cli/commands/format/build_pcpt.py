# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt format build-pcpt
==========================

Build a PCPT (targets / points of interest, ``.pcpt.json``) document
from a CSV or XLSX table -- see :class:`pycsamt.format.pointset.PointSet`.
"""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.options import no_color_option, overwrite_option, verbose_option
from ._base import _rich_table, fmt


@fmt.command("build-pcpt")
@click.argument(
    "source",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    metavar="SOURCE",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    metavar="OUT.pcpt.json",
    help="Output path (default: SOURCE with a .pcpt.json extension).",
)
@click.option(
    "--sheet", default=None, help="XLSX only: sheet name or index."
)
@click.option(
    "--header-row", type=int, default=None, help="XLSX only: header row index."
)
@click.option("--crs", default=None, help="CRS of the point coordinates.")
@click.option(
    "--document-id", default=None, help="Explicit document id (default: auto)."
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
@verbose_option
@no_color_option
@overwrite_option
def build_pcpt(
    source: Path,
    output_path: Path | None,
    sheet: str | None,
    header_row: int | None,
    crs: str | None,
    document_id: str | None,
    output_format: str,
    verbose: int,
    no_color: bool,
    overwrite: bool,
) -> None:
    """Build a PCPT points-of-interest document from SOURCE (.csv/.xlsx).

    \b
    Examples:
      pycsamt format build-pcpt targets.csv
      pycsamt format build-pcpt targets.xlsx --sheet Targets -o t.pcpt.json
    """
    from pycsamt.format.pointset import (
        points_from_csv,
        points_from_xlsx,
        write_points,
    )

    dst = output_path or source.with_suffix("").with_suffix(".pcpt.json")
    if dst.exists() and not overwrite:
        raise click.ClickException(
            f"{dst} exists — pass --overwrite to replace it."
        )

    suffix = source.suffix.lower()
    try:
        if suffix in (".xlsx", ".xlsm"):
            sheet_arg: str | int | None = sheet
            if sheet is not None and sheet.isdigit():
                sheet_arg = int(sheet)
            point_set = points_from_xlsx(
                source,
                sheet=sheet_arg,
                header_row=header_row,
                crs=crs,
                document_id=document_id,
            )
        else:
            point_set = points_from_csv(
                source, crs=crs, document_id=document_id
            )
    except Exception as exc:  # noqa: BLE001
        raise click.ClickException(f"Failed to build PCPT: {exc}") from exc

    written = write_points(point_set, dst)
    report = {
        "file": str(written),
        "document_id": point_set.document_id,
        "n_points": len(point_set.points),
        "crs": point_set.crs,
    }

    if output_format == "json":
        click.echo(json.dumps(report, indent=2, default=str))
        return
    _rich_table(
        "format build-pcpt",
        [
            ("wrote", report["file"]),
            ("document_id", report["document_id"]),
            ("points", str(report["n_points"])),
            ("crs", report["crs"] or "-"),
        ],
    )
    click.echo(f"\n✓ {written}")
