# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt pipe plugins — discover and list third-party pipeline steps."""

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


@pipe.command("plugins")
@click.option(
    "--group",
    default=None,
    metavar="ENTRY_POINT_GROUP",
    help="Override the entry-point group to scan (default: pycsamt.pipeline.steps).",
)
@click.option(
    "--strict",
    is_flag=True,
    default=False,
    help="Exit with a non-zero status if any plugin failed to load.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def plugins(
    ctx: click.Context,
    group: str | None,
    strict: bool,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Discover and list third-party pipeline steps.

    A plugin package registers pipeline steps by declaring a callable under
    the ``pycsamt.pipeline.steps`` entry-point group in its own
    ``pyproject.toml``.  This command loads every such entry point (calling
    each plugin's registration function) and reports the outcome, then
    lists every step currently registered with ``origin=plugin``.

    This is the one ``pycsamt pipe`` command that always discovers plugins
    — that is its job.  Other subcommands never do this implicitly (see
    ``pycsamt pipe --help`` for the ``--with-plugins`` opt-in flag), since
    scanning installed packages for entry points can take several seconds
    in a large environment.

    \b
    Examples:
      # List discovered plugins and what they registered
      pycsamt pipe plugins

      # Exit non-zero if any plugin failed to load (e.g. in CI)
      pycsamt pipe plugins --strict

      # Machine-readable
      pycsamt pipe plugins --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from pycsamt.pipeline import (  # noqa: PLC0415
        ENTRY_POINT_GROUP,
        STEP_REGISTRY,
        discover_plugins,
    )

    results = discover_plugins(group=group or ENTRY_POINT_GROUP, on_error="warn")
    plugin_steps = [s for s in STEP_REGISTRY.values() if s.origin == "plugin"]

    if output_format == "json":
        data = {
            "plugins": [
                {"name": r.name, "ok": r.ok, "error": r.error} for r in results
            ],
            "steps": [
                {
                    "code": s.code,
                    "name": s.name,
                    "label": s.label,
                    "category": s.category,
                }
                for s in plugin_steps
            ],
        }
        click.echo(json.dumps(data, indent=2))

    elif output_format == "csv":
        click.echo("plugin_name,ok,error")
        for r in results:
            err = (r.error or "").replace(",", ";")
            click.echo(f"{r.name},{r.ok},{err}")
        click.echo()
        click.echo("code,name,label,category")
        for s in plugin_steps:
            click.echo(f"{s.code},{s.name},{s.label},{s.category}")

    else:
        # text
        if not results:
            click.echo("No pipeline plugins found.")
        else:
            click.echo(f"Discovered {len(results)} pipeline plugin(s):")
            for r in results:
                status = "ok" if r.ok else f"FAILED ({r.error})"
                click.echo(f"  {r.name:<24} {status}")

        click.echo()
        if not plugin_steps:
            click.echo("No plugin steps registered.")
        else:
            click.echo(f"Registered {len(plugin_steps)} plugin step(s):")
            for s in plugin_steps:
                click.echo(f"  {s.code:<10} {s.name:<28} [{s.category}]  {s.label}")

    failed = [r.name for r in results if not r.ok]
    if strict and failed:
        raise click.ClickException(
            f"{len(failed)} plugin(s) failed to load: {', '.join(failed)}"
        )
