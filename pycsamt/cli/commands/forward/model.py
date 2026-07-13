# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt forward model — inspect and list geological earth model scenarios.

Sub-commands
------------
geology     List or inspect pre-defined geological scenarios from the
            pyCSAMT metadata catalog (``GEOLOGY_PRIORS``).

Usage
-----
::

    # list all scenarios (table view)
    pycsamt forward model geology

    # inspect one scenario
    pycsamt forward model geology --name sedimentary

    # JSON output for scripting
    pycsamt forward model geology --format json
    pycsamt forward model geology --name geothermal --format json
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
from ._base import forward

# ---------------------------------------------------------------------------
# model — sub-group
# ---------------------------------------------------------------------------


@forward.group("model")
@click.pass_context
def model(ctx: click.Context) -> None:
    """Inspect earth model scenarios used by the forward solvers.

    \b
    Sub-commands:
      geology     List or show pre-defined geological scenarios

    \b
    Examples:
      pycsamt forward model geology
      pycsamt forward model geology --name sedimentary
      pycsamt forward model geology --format json
    """
    ctx.ensure_object(dict)


# ---------------------------------------------------------------------------
# geology
# ---------------------------------------------------------------------------


@model.command("geology")
@click.option(
    "--name",
    "-n",
    default=None,
    metavar="SCENARIO",
    help=(
        "Name of the scenario to show in detail.  "
        "When omitted all available scenarios are listed."
    ),
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def model_geology(
    ctx: click.Context,
    name: str | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """List or inspect pre-defined geological model scenarios.

    Without ``--name`` prints a summary table of all built-in scenarios.
    With ``--name`` prints the full parameter details for that scenario,
    including the resistivity range, investigation depth, and n_layers hint.

    \b
    Built-in scenarios:
      sedimentary   crystalline   geothermal   marine   permafrost

    \b
    Examples:
      pycsamt forward model geology
      pycsamt forward model geology --name sedimentary
      pycsamt forward model geology --name geothermal --format json
      pycsamt forward model geology --format csv > scenarios.csv
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from ....forward.synthetic import (
        GEOLOGY_PRIORS,  # noqa: PLC0415
    )

    try:
        if name is not None:
            key = name.lower().strip()
            if key not in GEOLOGY_PRIORS:
                available = ", ".join(sorted(GEOLOGY_PRIORS))
                raise click.BadParameter(
                    f"{name!r} is not a known scenario.  "
                    f"Available: {available}",
                    param_hint="--name",
                )
            _emit_detail(key, GEOLOGY_PRIORS[key], output_format)
        else:
            _emit_table(GEOLOGY_PRIORS, output_format)
    except (click.exceptions.BadParameter, click.exceptions.UsageError):
        raise
    except Exception as exc:  # noqa: BLE001
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def _emit_table(priors: dict, output_format: str) -> None:
    import pandas as pd  # noqa: PLC0415

    rows = []
    for scenario, p in sorted(priors.items()):
        n_lo, n_hi = p.get("n_layers", ("—", "—"))
        lr_lo, lr_hi = p.get("log_rho_range", ("—", "—"))
        d_lo, d_hi = p.get("depth_max_range", ("—", "—"))
        rows.append(
            {
                "scenario": scenario,
                "n_layers": f"{n_lo}–{n_hi}",
                "log10_rho_range": f"{lr_lo:.1f}–{lr_hi:.1f}",
                "depth_max_m": f"{d_lo:.0f}–{d_hi:.0f}",
                "description": (p.get("description", "") or "")[:60],
            }
        )

    df = pd.DataFrame(rows)

    if output_format == "json":
        click.echo(
            df.to_json(orient="records", indent=2, default_handler=str)
        )
    elif output_format == "csv":
        click.echo(df.to_csv(index=False))
    else:
        click.echo(df.to_string(index=False))


def _emit_detail(name: str, prior: dict, output_format: str) -> None:
    if output_format == "json":
        click.echo(
            json.dumps({"scenario": name, **prior}, indent=2, default=str)
        )
        return

    if output_format == "csv":
        import pandas as pd  # noqa: PLC0415

        rows = [{"key": k, "value": str(v)} for k, v in prior.items()]
        click.echo(pd.DataFrame(rows).to_csv(index=False))
        return

    n_lo, n_hi = prior.get("n_layers", ("—", "—"))
    lr_lo, lr_hi = prior.get("log_rho_range", ("—", "—"))
    d_lo, d_hi = prior.get("depth_max_range", ("—", "—"))
    desc = prior.get("description", "")

    click.echo(f"\nScenario      : {name}")
    click.echo(f"n_layers      : {n_lo} – {n_hi}")
    click.echo(
        f"log10(ρ) range: {lr_lo:.2f} – {lr_hi:.2f}  "
        f"(ρ ≈ {10**lr_lo:.0f} – {10**lr_hi:.0f} Ω·m)"
    )
    click.echo(f"depth_max     : {d_lo:.0f} – {d_hi:.0f} m")
    if desc:
        click.echo(f"description   : {desc}")
