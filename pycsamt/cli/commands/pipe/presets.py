# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt pipe presets — list and expand named pipeline presets."""

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


@pipe.command("presets")
@click.option(
    "--expand",
    "preset_name",
    default=None,
    metavar="NAME",
    help=(
        "Show the full step list for a specific preset.  "
        "Example: --expand full_processing"
    ),
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def presets(
    ctx: click.Context,
    preset_name: str | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """List available pipeline presets.

    Presets are named, opinionated processing workflows that cover
    the most common AMT/MT use cases.  Any preset can be used as
    the starting point for a run::

        pycsamt pipe run --preset full_processing

    or exported to a config file for further customisation::

        pycsamt pipe init --preset full_processing --format yaml

    \b
    Examples:
      # All presets with descriptions
      pycsamt pipe presets

      # Expand a specific preset: show its steps
      pycsamt pipe presets --expand full_processing

      # Machine-readable
      pycsamt pipe presets --format json
      pycsamt pipe presets --expand basic_qc --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from pycsamt.pipeline import (  # noqa: PLC0415
        Pipeline,
        get_preset,
        list_presets,
    )

    # ── Expand one preset ─────────────────────────────────────────────────
    if preset_name is not None:
        try:
            preset = get_preset(preset_name)
        except KeyError:
            available = [p.name for p in list_presets()]
            raise click.BadParameter(
                f"Unknown preset {preset_name!r}.  Available: {available}",
                param_hint="--expand",
            )

        pipe_obj = Pipeline.from_preset(preset_name)

        if output_format == "json":
            data = {
                "name": preset.name,
                "description": preset.description,
                "n_steps": len(preset.steps),
                "steps": [
                    {
                        "label": lbl,
                        "code": step.spec.code,
                        "name": step.spec.name,
                        "category": step.spec.category,
                        "params": {
                            k: list(v) if isinstance(v, tuple) else v
                            for k, v in step.params.items()
                        },
                    }
                    for lbl, step in preset.steps
                ],
            }
            click.echo(json.dumps(data, indent=2))
        else:
            click.echo(str(pipe_obj))
            click.echo()
            click.echo(f"  {preset.description}")
        return

    # ── List all presets ──────────────────────────────────────────────────
    all_presets = list_presets()

    if output_format == "json":
        data = [
            {
                "name": p.name,
                "description": p.description,
                "n_steps": len(p.steps),
                "codes": [s.spec.code for _, s in p.steps],
            }
            for p in all_presets
        ]
        click.echo(json.dumps(data, indent=2))

    elif output_format == "csv":
        click.echo("name,description,n_steps,codes")
        for p in all_presets:
            codes = "|".join(s.spec.code for _, s in p.steps)
            click.echo(f"{p.name},{p.description!r},{len(p.steps)},{codes}")

    else:
        from pycsamt.pipeline import (
            preset_catalogue,  # noqa: PLC0415
        )

        click.echo(preset_catalogue())
        click.echo(
            "  Tip: pycsamt pipe presets --expand <name>  to see a preset's steps."
        )
        click.echo("       pycsamt pipe run --preset <name>       to run a preset.")
