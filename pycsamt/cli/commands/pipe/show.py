# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt pipe show — pretty-print a pipeline config file or preset."""

from __future__ import annotations

import json
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    no_color_option,
    verbose_option,
)
from ._base import pipe


@pipe.command("show")
@click.argument(
    "config_file",
    type=click.Path(exists=True, path_type=Path),
    metavar="[CONFIG_FILE]",
    required=False,
    default=None,
)
@click.option(
    "--preset",
    default=None,
    metavar="NAME",
    help="Show a named preset instead of a file.",
)
@click.option(
    "--from-step",
    default=None,
    metavar="LABEL_OR_CODE",
    help="Preview the pipeline after slicing from this step.",
)
@click.option(
    "--until-step",
    default=None,
    metavar="LABEL_OR_CODE",
    help="Preview the pipeline after slicing until this step.",
)
@click.option(
    "--n-steps",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Preview only the first N steps.",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def show(
    ctx: click.Context,
    config_file: Path | None,
    preset: str | None,
    from_step: str | None,
    until_step: str | None,
    n_steps: int | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Pretty-print a pipeline config file or named preset.

    Loads the pipeline and renders it as the same formatted table
    produced by ``print(pipeline)``.  Combine with --from-step /
    --until-step / --n-steps to preview a sliced version before
    passing those flags to  pycsamt pipe run.

    \b
    Examples:
      # Inspect a YAML config
      pycsamt pipe show workflow.yaml

      # Inspect a preset
      pycsamt pipe show --preset publication_ready

      # Preview what the first 3 steps look like
      pycsamt pipe show workflow.yaml --n-steps 3

      # Machine-readable step list
      pycsamt pipe show --preset full_processing --format json
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from ._base import (  # noqa: PLC0415
        _resolve_pipeline,
        _trim_pipeline,
    )

    if config_file is None and preset is None:
        raise click.UsageError(
            "Provide a CONFIG_FILE argument or --preset NAME.\n"
            "Run  pycsamt pipe presets  to list preset names."
        )

    try:
        pipeline = _resolve_pipeline(
            config=config_file,
            preset=preset,
            steps=None,
            name=None,
            verbose=verbose,
        )
    except Exception as exc:
        raise click.ClickException(f"Failed to load: {exc}") from exc

    pipeline = _trim_pipeline(pipeline, from_step, until_step, n_steps)

    if output_format == "json":
        data = {
            "name": pipeline.name,
            "n_steps": len(pipeline),
            "steps": [
                {
                    "label": lbl,
                    "code": step.spec.code,
                    "name": step.spec.name,
                    "category": step.spec.category,
                    "label_long": step.spec.label,
                    "params": {
                        k: list(v) if isinstance(v, tuple) else v
                        for k, v in step.params.items()
                    },
                }
                for lbl, step in pipeline._steps
            ],
        }
        click.echo(json.dumps(data, indent=2))

    elif output_format == "csv":
        click.echo("idx,label,code,name,category,params")
        for i, (lbl, step) in enumerate(pipeline._steps, start=1):
            params_str = ";".join(
                f"{k}={v}" for k, v in step.params.items()
            )
            click.echo(
                f"{i},{lbl},{step.spec.code},"
                f"{step.spec.name},{step.spec.category},{params_str}"
            )

    else:
        click.echo(str(pipeline))
