# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt jones validate — structural validation of Jones J-files."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    no_color_option,
    verbose_option,
)
from ._base import jones


@jones.command("validate")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="J_OR_DIR",
)
@click.option(
    "--deep/--no-deep",
    default=True,
    show_default=True,
    help=(
        "Deep validation scans the file content for valid data blocks "
        "(station + datatype + count + data rows).  "
        "--no-deep checks the file extension only."
    ),
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
@verbose_option
@no_color_option
@click.pass_context
def validate(
    ctx: click.Context,
    source: Path,
    deep: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Check Jones J-files for structural validity.

    Each file is tested with :func:`~pycsamt.jones.validation.is_j_file`.
    Deep mode (default) scans the content for at least one complete data
    block (station + datatype + count + rows).  ``--no-deep`` only checks
    the file extension.

    Exit code is 0 when all files pass, 1 when any fail.

    \b
    Examples:
      pycsamt jones validate data/j/
      pycsamt jones validate data/j/S01.j
      pycsamt jones validate data/j/ --no-deep
      pycsamt jones validate data/j/ --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from pycsamt.jones.validation import (
        is_j_file,  # noqa: PLC0415
    )

    # Collect target files — accept any extension used for J-files
    _J_EXTS = {".j", ".dat", ".txt"}

    if source.is_file():
        files = [source]
    else:
        files = sorted(
            f
            for f in source.rglob("*")
            if f.suffix.lower() in _J_EXTS and f.is_file()
        )

    if not files:
        click.echo(
            f"No J-files (.j / .dat / .txt) found in {source}.", err=True
        )
        sys.exit(1)

    results: list[dict] = []
    for f in files:
        try:
            is_j_file(f, deep=deep)
            results.append({"path": str(f), "valid": True, "error": None})
        except Exception as exc:  # noqa: BLE001
            results.append({"path": str(f), "valid": False, "error": str(exc)})

    n_ok = sum(1 for r in results if r["valid"])
    n_fail = len(results) - n_ok
    mode = "deep" if deep else "extension-only"

    if output_format == "json":
        click.echo(
            json.dumps(
                {
                    "n_files": len(results),
                    "n_ok": n_ok,
                    "n_fail": n_fail,
                    "mode": mode,
                    "results": results,
                },
                indent=2,
            )
        )
    else:
        click.echo(
            f"Validated {len(results)} file(s)  [{mode}]  "
            f"→  {n_ok} ok  /  {n_fail} failed"
        )
        if n_fail:
            click.echo()
            for r in results:
                if not r["valid"]:
                    click.echo(f"  ✗  {r['path']}")
                    if r["error"]:
                        click.echo(
                            f"       {r['error'].splitlines()[0]}", err=True
                        )
        elif len(results) <= 10:
            for r in results:
                click.echo(f"  ✓  {r['path']}")

    if n_fail > 0:
        sys.exit(1)
