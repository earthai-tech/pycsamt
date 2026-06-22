# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt edi validate — structural validation of EDI files."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import no_color_option, verbose_option
from ._base import edi


@edi.command("validate")
@click.argument(
    "source",
    type=click.Path(exists=True, path_type=Path),
    metavar="EDI_OR_DIR",
)
@click.option(
    "--deep/--no-deep",
    default=True,
    show_default=True,
    help=(
        "Deep validation scans file content for structural tags "
        "(>HEAD, >FREQ, >=MTSECT, >END etc.).  "
        "--no-deep checks the .edi extension only."
    ),
)
@click.option(
    "--format", "output_format",
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
    """Check EDI files for structural validity.

    Each file is tested with :class:`~pycsamt.seg.validation.IsEdi`.
    The ``--deep`` mode (default) scans the file content for required
    structural tags; ``--no-deep`` only checks the ``.edi`` extension.

    Exit code is 0 when all files pass, 1 when any fail.

    \b
    Examples:
      pycsamt edi validate data/3edis/
      pycsamt edi validate data/AMT/ --deep
      pycsamt edi validate data/3edis/HBH03.edi
      pycsamt edi validate data/AMT/ --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from pycsamt.seg.validation import IsEdi  # noqa: PLC0415

    # Collect target files
    if source.is_file():
        files = [source]
    else:
        files = sorted(source.rglob("*.edi"))

    if not files:
        click.echo(f"No .edi files found in {source}.", err=True)
        sys.exit(1)

    results: list[dict] = []
    for f in files:
        try:
            IsEdi._assert_edi(f, deep=deep)
            results.append({"path": str(f), "valid": True, "error": None})
        except Exception as exc:  # noqa: BLE001
            results.append({"path": str(f), "valid": False, "error": str(exc)})

    n_ok   = sum(1 for r in results if r["valid"])
    n_fail = len(results) - n_ok

    if output_format == "json":
        click.echo(json.dumps({
            "n_files": len(results),
            "n_ok":    n_ok,
            "n_fail":  n_fail,
            "results": results,
        }, indent=2))
    else:
        _print_text(results, n_ok, n_fail, deep)

    if n_fail > 0:
        sys.exit(1)


def _print_text(results: list[dict], n_ok: int, n_fail: int, deep: bool) -> None:
    mode = "deep" if deep else "extension-only"
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
                    first_line = r["error"].splitlines()[0]
                    click.echo(f"       {first_line}", err=True)
    elif len(results) <= 10:
        for r in results:
            click.echo(f"  ✓  {r['path']}")
