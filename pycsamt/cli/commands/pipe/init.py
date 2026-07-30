# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt pipe init — scaffold a ready-to-edit pipeline config file."""

from __future__ import annotations

from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    no_color_option,
    verbose_option,
)
from ._base import pipe


@pipe.command("init")
@click.option(
    "--preset",
    default=None,
    metavar="NAME",
    help=(
        "Seed the scaffold with steps from this preset.  "
        "All other available steps appear as comments.  "
        "Run  pycsamt pipe presets  to list names."
    ),
)
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["yaml", "json", "py"], case_sensitive=False),
    default="yaml",
    show_default=True,
    help="Output config format.",
)
@click.option(
    "--name",
    default="my_workflow",
    show_default=True,
    metavar="TEXT",
    help="Pipeline name written into the config.",
)
@click.option(
    "--outdir",
    "pipe_outdir",
    default="pipe_results",
    show_default=True,
    metavar="DIR",
    help="Default output directory embedded in the config.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    metavar="FILE_OR_DIR",
    help=(
        "Path to write the generated config.  "
        "If a directory is given the filename is derived from --name and --format.  "
        "Defaults to <name>.<format> in the current directory."
    ),
)
@click.option(
    "--print",
    "print_only",
    is_flag=True,
    default=False,
    help="Print the generated config to stdout instead of writing to disk.",
)
@verbose_option
@no_color_option
@click.pass_context
def init(
    ctx: click.Context,
    preset: str | None,
    fmt: str,
    name: str,
    pipe_outdir: str,
    output_path: Path | None,
    print_only: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Generate a ready-to-edit pipeline configuration file.

    The scaffold contains the active steps (from --preset or a sensible
    default set) uncommented, with all other available steps shown as
    comments so you can easily add them.

    \b
    Format options:
      yaml  (default)  — human-readable, most common
      json             — machine-friendly
      py               — Python dict, most editable; supports inline logic

    \b
    Examples:
      # Default YAML scaffold in current directory
      pycsamt pipe init

      # YAML from a preset
      pycsamt pipe init --preset full_processing --name willy_survey

      # Python format, custom output path
      pycsamt pipe init --format py --preset basic_qc -o config/willy.py

      # JSON format, print to stdout for inspection
      pycsamt pipe init --format json --print

      # Custom name and output dir embedded in config
      pycsamt pipe init --name l22_profile --outdir /data/results/ -o config/
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from pycsamt.pipeline import Pipeline  # noqa: PLC0415

    # Resolve output path
    if not print_only:
        if output_path is None:
            output_path = Path(f"{name}.{fmt.lower()}")
        elif output_path.is_dir() or (
            not output_path.suffix and not output_path.exists()
        ):
            output_path = output_path / f"{name}.{fmt.lower()}"

    try:
        content = Pipeline.scaffold(
            path=None if print_only else output_path,
            fmt=fmt.lower(),
            preset=preset,
            name=name,
            outdir=pipe_outdir,
        )
    except KeyError as exc:
        raise click.ClickException(str(exc)) from exc
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="--format") from exc

    if print_only:
        click.echo(content)
    else:
        click.echo(f"Created: {output_path}")
        if verbose >= 1:
            click.echo()
            click.echo("  Load with:")
            if fmt.lower() == "yaml":
                click.echo(
                    f"    pipe = Pipeline.from_yaml({str(output_path)!r})"
                )
            elif fmt.lower() == "json":
                click.echo(
                    f"    pipe = Pipeline.from_json({str(output_path)!r})"
                )
            else:
                click.echo(
                    f"    pipe = Pipeline.from_py({str(output_path)!r})"
                )
            click.echo()
            click.echo("  Run with:")
            click.echo(f"    pycsamt pipe run --config {output_path}")
