# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt interp rocks — browse or query the built-in rock-physics database.

Without ``--rho`` prints the full 25-entry database sorted by resistivity.
With ``--rho VALUE`` classifies that linear Ω·m value and returns the
best-matching geological unit.

Usage
-----
::

    # list all rocks
    pycsamt interp rocks

    # classify a specific resistivity value
    pycsamt interp rocks --rho 250

    # JSON output for scripting
    pycsamt interp rocks --format json
    pycsamt interp rocks --rho 250 --format json
"""

from __future__ import annotations

import json
import sys

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    no_color_option,
    verbose_option,
)
from ._base import interp


@interp.command("rocks")
@click.option(
    "--rho",
    "-r",
    type=float,
    default=None,
    metavar="OHM_M",
    help=(
        "Resistivity value to classify in Ω·m (linear).  "
        "Omit to list the full database."
    ),
)
@click.option(
    "--method",
    type=click.Choice(["nearest", "overlap"], case_sensitive=False),
    default="nearest",
    show_default=True,
    help=(
        "Classification method.  "
        "'nearest': closest geometric-mean resistivity.  "
        "'overlap': first entry whose range brackets the value."
    ),
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def rocks(
    ctx: click.Context,
    rho: float | None,
    method: str,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Browse or query the built-in rock-physics database.

    Without ``--rho``, prints a sorted table of all 25 rock entries with
    their resistivity ranges and descriptions.  With ``--rho``, returns
    the best-matching geological unit for that resistivity value.

    \b
    Examples:
      pycsamt interp rocks
      pycsamt interp rocks --rho 250
      pycsamt interp rocks --rho 0.05 --method overlap
      pycsamt interp rocks --format json
      pycsamt interp rocks --rho 5000 --format csv
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from ....interp.lithology import (
        RockDatabase,  # noqa: PLC0415
    )

    try:
        db = RockDatabase.default()

        if rho is not None:
            if rho <= 0:
                raise click.BadParameter(
                    "Resistivity must be a positive value (Ω·m).",
                    param_hint="--rho",
                )
            entry = db.classify(rho, method=method)
            _emit_entry(rho, entry, output_format)
        else:
            _emit_database(db, output_format)

    except (click.exceptions.BadParameter, click.exceptions.UsageError):
        raise
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def _emit_database(db, output_format: str) -> None:
    import pandas as pd  # noqa: PLC0415

    rows = [
        {
            "code": e.code,
            "name": e.name,
            "rho_min": e.rho_min,
            "rho_max": e.rho_max,
            "rho_mid": round(e.rho_mid, 2),
            "description": e.description,
        }
        for e in db._entries
    ]
    df = pd.DataFrame(rows).sort_values("rho_min").reset_index(drop=True)

    if output_format == "json":
        click.echo(df.to_json(orient="records", indent=2, default_handler=str))
    elif output_format == "csv":
        click.echo(df.to_csv(index=False))
    else:
        click.echo(df.to_string(index=False))


def _emit_entry(rho: float, entry, output_format: str) -> None:
    data = {
        "query_rho_Ohm_m": rho,
        "code": entry.code,
        "name": entry.name,
        "rho_min_Ohm_m": entry.rho_min,
        "rho_max_Ohm_m": entry.rho_max,
        "rho_mid_Ohm_m": round(entry.rho_mid, 2),
        "description": entry.description,
    }
    if output_format == "json":
        click.echo(json.dumps(data, indent=2, default=str))
    elif output_format == "csv":
        import pandas as pd  # noqa: PLC0415

        click.echo(pd.DataFrame([data]).to_csv(index=False))
    else:
        click.echo(f"\nQuery  : {rho} Ω·m")
        click.echo(f"Match  : {entry.name}  (code {entry.code})")
        click.echo(
            f"Range  : {entry.rho_min} – {entry.rho_max} Ω·m  "
            f"(mid ≈ {entry.rho_mid:.1f} Ω·m)"
        )
        if entry.description:
            click.echo(f"Note   : {entry.description}")
