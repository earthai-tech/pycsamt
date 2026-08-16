# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt pipe steps — browse and inspect available processing steps."""

from __future__ import annotations

import json

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    no_color_option,
    verbose_option,
)
from ._base import pipe


@pipe.command("steps")
@click.option(
    "--category",
    default=None,
    metavar="NAME",
    help=(
        "Filter by step category.  "
        "One of: frequency, noise_removal, static_shift, tensor, "
        "dimensionality, skew, source_effects, qc, preview, export."
    ),
)
@click.option(
    "--info",
    "step_code",
    default=None,
    metavar="CODE_OR_NAME",
    help=(
        "Show detailed information for a single step.  "
        "Accepts both code (NR001) and snake-case name (notch_powerline)."
    ),
)
@click.option(
    "--codes-only",
    is_flag=True,
    default=False,
    help="Print only step codes, one per line (useful for scripting).",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def steps(
    ctx: click.Context,
    category: str | None,
    step_code: str | None,
    codes_only: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """List and inspect available MT processing steps.

    Without options, prints a grouped catalogue of all 55 registered
    steps.  Use --category to narrow the view, or --info CODE to see
    full details about a specific step.

    \b
    Examples:
      # All steps, grouped by category
      pycsamt pipe steps

      # Only noise-removal steps
      pycsamt pipe steps --category noise_removal

      # Full info for one step
      pycsamt pipe steps --info NR001
      pycsamt pipe steps --info notch_powerline

      # Machine-readable list of all codes
      pycsamt pipe steps --codes-only --format json

      # Pipe codes into another command
      pycsamt pipe steps --category static_shift --codes-only
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from pycsamt.pipeline import (  # noqa: PLC0415
        Pipeline,
        categories,
        list_steps,
        lookup_step,
    )

    # ── --info: single step detail ────────────────────────────────────────
    if step_code is not None:
        try:
            spec = lookup_step(step_code)
        except KeyError:
            raise click.BadParameter(
                f"Unknown step {step_code!r}.  "
                "Use  pycsamt pipe steps  to list all codes.",
                param_hint="--info",
            )

        if output_format == "json":
            qc_names = [fn for _, fn in spec.qc_defs]
            data = {
                "code": spec.code,
                "name": spec.name,
                "label": spec.label,
                "category": spec.category,
                "module": spec.mod,
                "function": spec.fn_name,
                "defaults": {
                    k: list(v) if isinstance(v, tuple) else v
                    for k, v in spec.defaults.items()
                },
                "qc_plots": qc_names,
                "returns_sites": spec.returns_sites,
            }
            click.echo(json.dumps(data, indent=2))
        else:
            click.echo(Pipeline.step_info(step_code))
        return

    # ── --codes-only ──────────────────────────────────────────────────────
    if codes_only:
        specs = list_steps(category)
        if output_format == "json":
            click.echo(json.dumps([s.code for s in specs]))
        else:
            for s in specs:
                click.echo(s.code)
        return

    # ── Full catalogue ────────────────────────────────────────────────────
    if output_format == "json":
        cats = [category] if category else categories()
        data = {}
        for cat in cats:
            specs = list_steps(cat)
            data[cat] = [
                {
                    "code": s.code,
                    "name": s.name,
                    "label": s.label,
                    "defaults": {
                        k: list(v) if isinstance(v, tuple) else v
                        for k, v in s.defaults.items()
                    },
                    "returns_sites": s.returns_sites,
                }
                for s in specs
            ]
        click.echo(json.dumps(data, indent=2))
    elif output_format == "csv":
        click.echo("code,name,label,category,returns_sites")
        for s in list_steps(category):
            click.echo(
                f"{s.code},{s.name},{s.label},{s.category},{s.returns_sites}"
            )
    else:
        click.echo(Pipeline.catalogue(category))
