# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt format detect
=====================

Report what a path is — a PCSF/PCSM file, a classical-solver working
directory, or an AI array bundle — and how ``pycsamt format convert``
would treat it, without touching heavy array data.
"""

from __future__ import annotations

import json
from pathlib import Path

import click

from ._base import _detect, _rich_table, fmt


@fmt.command("detect")
@click.argument(
    "path",
    type=click.Path(exists=True, path_type=Path),
    metavar="PATH",
)
@click.option(
    "--solver",
    type=click.Choice(["occam2d", "modem", "mare2dem"], case_sensitive=False),
    default=None,
    help="Force the solver backend instead of fingerprinting it.",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
def detect(path: Path, solver: str | None, output_format: str) -> None:
    """Classify PATH (a file or a solver working directory).

    \b
    Examples:
      pycsamt format detect data/occam2D
      pycsamt format detect data/mare2dem/demo_mt_inversion
      pycsamt format detect unet_prediction.npz -f json
    """
    sk = _detect(path, solver)

    if output_format == "json":
        click.echo(json.dumps(sk.to_dict(), indent=2, default=str))
        return

    rows = [
        ("path", str(sk.path)),
        ("category", sk.category),
        ("backend", sk.backend or "-"),
        ("geometry", sk.geometry or "-"),
        ("→ geometry", sk.target_geometry or "-"),
        ("confidence", sk.confidence),
        ("convertible", "yes" if sk.convertible else "no"),
        ("detail", sk.detail),
    ]
    for key, value in sk.hints.items():
        if value is None or key in {"keys", "shapes"}:
            continue
        rows.append((f"hint: {key}", str(value)))
    _rich_table("format detect", rows)

    if not sk.convertible:
        raise SystemExit(1)
