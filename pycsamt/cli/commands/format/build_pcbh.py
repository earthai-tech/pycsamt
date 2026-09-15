# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt format build-pcbh
==========================

Build a PCBH (borehole, ``.pcbh.json``) document from a combined CSV, a
directory of relational tables, an XLSX workbook, or a single LAS 2.0
log -- see :mod:`pycsamt.format.borehole`.
"""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.options import no_color_option, overwrite_option, verbose_option
from ._base import _rich_table, fmt

_FROM_CHOICES = ["auto", "csv", "csv-dir", "xlsx", "las"]


def _infer_from(source: Path) -> str:
    if source.is_dir():
        return "csv-dir"
    suffix = source.suffix.lower()
    if suffix in (".xlsx", ".xlsm"):
        return "xlsx"
    if suffix == ".las":
        return "las"
    return "csv"


@fmt.command("build-pcbh")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="SOURCE",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    metavar="OUT.pcbh.json",
    help=(
        "Output path (default: SOURCE with a .pcbh.json extension, or "
        "SOURCE/borehole.pcbh.json for a directory)."
    ),
)
@click.option(
    "--from",
    "source_kind",
    type=click.Choice(_FROM_CHOICES, case_sensitive=False),
    default="auto",
    show_default=True,
    help="Source kind (default: inferred from SOURCE's suffix / type).",
)
@click.option(
    "--collar-id", default=None, help="LAS only: this borehole's id."
)
@click.option("--x", type=float, default=None, help="LAS only: collar X.")
@click.option("--y", type=float, default=None, help="LAS only: collar Y.")
@click.option(
    "--z", type=float, default=None, help="LAS only: collar Z (elevation)."
)
@click.option(
    "--crs",
    "crs_horizontal",
    default=None,
    help="LAS only: horizontal CRS (e.g. 'EPSG:32650').",
)
@click.option(
    "--document-id", default=None, help="Explicit document id (default: auto)."
)
@click.option(
    "--created-by",
    default="pycsamt format build-pcbh",
    show_default=True,
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
def build_pcbh(
    source: Path,
    output_path: Path | None,
    source_kind: str,
    collar_id: str | None,
    x: float | None,
    y: float | None,
    z: float | None,
    crs_horizontal: str | None,
    document_id: str | None,
    created_by: str,
    output_format: str,
    verbose: int,
    no_color: bool,
    overwrite: bool,
) -> None:
    """Build a PCBH borehole document from SOURCE.

    SOURCE is a combined collar+interval CSV, a directory of relational
    tables (with an ``import.yaml`` manifest), an XLSX workbook, or a
    single LAS 2.0 log (needs --collar-id/--x/--y/--z/--crs).

    \b
    Examples:
      pycsamt format build-pcbh boreholes.csv
      pycsamt format build-pcbh project_dir/ --from csv-dir
      pycsamt format build-pcbh logs.xlsx -o site.pcbh.json
      pycsamt format build-pcbh ZK01.las --collar-id ZK01 \\
          --x 512340 --y 2894210 --z 118.4 --crs EPSG:32650
    """
    kind = source_kind if source_kind != "auto" else _infer_from(source)

    default_name = (
        (source / "borehole.pcbh.json")
        if source.is_dir()
        else source.with_suffix("").with_suffix(".pcbh.json")
    )
    dst = output_path or default_name
    if dst.exists() and not overwrite:
        raise click.ClickException(
            f"{dst} exists — pass --overwrite to replace it."
        )

    report = None
    try:
        if kind == "csv":
            from pycsamt.format.borehole.csvio import boreholes_from_csv

            document, report = boreholes_from_csv(
                source, document_id=document_id, created_by=created_by
            )
        elif kind == "csv-dir":
            from pycsamt.format.borehole.relational import (
                boreholes_from_csv_directory,
            )

            document, report = boreholes_from_csv_directory(
                source, document_id=document_id, created_by=created_by
            )
        elif kind == "xlsx":
            from pycsamt.format.borehole.xlsxio import boreholes_from_xlsx

            document, report = boreholes_from_xlsx(
                source, document_id=document_id, created_by=created_by
            )
        elif kind == "las":
            if collar_id is None or x is None or y is None or z is None or (
                crs_horizontal is None
            ):
                raise click.UsageError(
                    "--from las needs --collar-id, --x, --y, --z and --crs."
                )
            from pycsamt.format.borehole.lasio import borehole_from_las
            from pycsamt.format.borehole.schema import Collar

            collar = Collar(x=x, y=y, z=z)
            document, report = borehole_from_las(
                source, collar=collar, crs_horizontal=crs_horizontal
            )
            document.document_id = document_id or document.document_id
            document.created_by = created_by
        else:  # pragma: no cover - guarded by click.Choice
            raise click.UsageError(f"Unknown --from kind: {kind!r}")
    except click.ClickException:
        raise
    except Exception as exc:  # noqa: BLE001
        raise click.ClickException(f"Failed to build PCBH: {exc}") from exc

    from pycsamt.format.borehole.jsonio import write_pcbh

    written = write_pcbh(document, dst)
    result = {
        "file": str(written),
        "document_id": document.document_id,
        "n_boreholes": len(document.boreholes),
    }
    if report is not None:
        result["rows_read"] = report.rows_read
        result["rows_accepted"] = report.rows_accepted
        result["rows_rejected"] = report.rows_rejected
        result["n_issues"] = len(report.issues)

    if output_format == "json":
        click.echo(json.dumps(result, indent=2, default=str))
        return
    rows = [
        ("wrote", result["file"]),
        ("document_id", result["document_id"]),
        ("boreholes", str(result["n_boreholes"])),
    ]
    if report is not None:
        rows += [
            ("rows read/accepted/rejected", (
                f"{report.rows_read}/{report.rows_accepted}/"
                f"{report.rows_rejected}"
            )),
            ("issues", str(len(report.issues))),
        ]
    _rich_table("format build-pcbh", rows)
    click.echo(f"\n✓ {written}")
    if report is not None and report.issues and verbose:
        for issue in report.issues[:20]:
            where = (
                f"row {issue.row}" if issue.row is not None else "-"
            )
            click.echo(f"  [{issue.severity}] {where}: {issue.message}")
