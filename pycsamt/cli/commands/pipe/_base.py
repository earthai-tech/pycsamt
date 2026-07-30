# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt pipe — root Click group and shared helpers.

Every sub-module in this package imports ``pipe`` from here and
decorates its commands with ``@pipe.command(...)`` so that
``__init__.py`` only needs to import the sub-modules to trigger
registrations.

Shared utilities
----------------
_resolve_pipeline   Build a Pipeline from --config / --preset / --steps.
_resolve_sites      Resolve any EDI source to a Sites object.
_trim_pipeline      Apply --from-step / --until-step / --n-steps slicing.
_format_run_result  Render a PipelineResult as text / JSON / CSV.
_rich_pipe_table    Pretty-print a key/value summary (rich or plain).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import click

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _resolve_pipeline(
    config: Path | None,
    preset: str | None,
    steps: list[str] | None,
    name: str | None,
    verbose: int = 0,
) -> Any:
    """Build a :class:`~pycsamt.pipeline.Pipeline`.

    Priority: *config* > *preset* > *steps* (explicit codes).
    Raises :class:`click.UsageError` when none of the three is supplied.
    """
    from pycsamt.pipeline import (  # noqa: PLC0415
        Pipeline,
        Step,
    )

    if config is not None:
        suffix = config.suffix.lower()
        if suffix in {".yaml", ".yml"}:
            pipe = Pipeline.from_yaml(config)
        elif suffix == ".json":
            pipe = Pipeline.from_json(config)
        elif suffix == ".py":
            pipe = Pipeline.from_py(config)
        else:
            raise click.BadParameter(
                f"Unsupported config format {suffix!r}.  "
                "Expected .yaml, .json, or .py",
                param_hint="--config",
            )
        if name:
            pipe.name = name
        if verbose >= 1:
            click.echo(f"  Loaded pipeline from {config}", err=True)
        return pipe

    if preset is not None:
        try:
            pipe = Pipeline.from_preset(preset)
        except KeyError:
            from pycsamt.pipeline import (
                list_presets,  # noqa: PLC0415
            )

            names = [p.name for p in list_presets()]
            raise click.BadParameter(
                f"Unknown preset {preset!r}.  Available: {names}",
                param_hint="--preset",
            )
        if name:
            pipe.name = name
        if verbose >= 1:
            click.echo(
                f"  Loaded preset {preset!r} ({len(pipe)} steps)", err=True
            )
        return pipe

    if steps:
        from pycsamt.pipeline import (
            lookup_step,  # noqa: PLC0415
        )

        step_objects = [(lookup_step(code).name, Step(code)) for code in steps]
        pipe = Pipeline(step_objects, name=name or "cli_pipeline")
        if verbose >= 1:
            click.echo(
                f"  Built pipeline from --steps: {', '.join(steps)}",
                err=True,
            )
        return pipe

    raise click.UsageError(
        "No pipeline specified.  Provide one of:\n"
        "  --config FILE    (YAML/JSON/Python config)\n"
        "  --preset NAME    (built-in preset)\n"
        "  --steps NR001,FREQ001,...  (explicit step codes)\n\n"
        "Run  pycsamt pipe presets  to list presets, or\n"
        "     pycsamt pipe steps    to list all step codes."
    )


def _resolve_sites(
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    verbose: int,
) -> Any:
    """Resolve any combination of positional path / --survey / active context."""
    from pycsamt.cli.survey import (
        resolve_survey,  # noqa: PLC0415
    )

    return resolve_survey(
        edi_source or survey_path,
        fresh=fresh,
        verbose=verbose,
    )


def _trim_pipeline(
    pipe: Any,
    from_step: str | None,
    until_step: str | None,
    n_steps: int | None,
) -> Any:
    """Slice *pipe._steps* according to --from-step / --until-step / --n-steps.

    All three are optional and composable:
    - ``--from-step align``  → skip steps before the one labelled/coded "align"
    - ``--until-step SS001`` → stop after the step with code SS001
    - ``--n-steps 3``        → keep only the first 3 steps
    """
    steps = list(pipe._steps)

    if from_step:
        idx = _find_step_index(steps, from_step)
        if idx is None:
            raise click.BadParameter(
                f"No step labelled or coded {from_step!r} in this pipeline.",
                param_hint="--from-step",
            )
        steps = steps[idx:]

    if until_step:
        idx = _find_step_index(steps, until_step)
        if idx is None:
            raise click.BadParameter(
                f"No step labelled or coded {until_step!r} in this pipeline.",
                param_hint="--until-step",
            )
        steps = steps[: idx + 1]

    if n_steps is not None and n_steps > 0:
        steps = steps[:n_steps]

    pipe._steps = steps
    return pipe


def _find_step_index(
    steps: list[tuple[str, Any]],
    query: str,
) -> int | None:
    """Return the 0-based index of the step matching *query* (label or code)."""
    for i, (label, step) in enumerate(steps):
        if (
            label == query
            or step.spec.code == query
            or step.spec.name == query
        ):
            return i
    return None


def _format_run_result(
    result: Any,
    output_format: str,
) -> str:
    """Render a :class:`~pycsamt.pipeline.PipelineResult` as text/JSON/CSV."""
    if output_format == "json":
        data: dict[str, Any] = {
            "pipeline": result.pipeline_name,
            "ok": result.ok,
            "n_errors": result.n_errors,
            "elapsed_sec": round(result.elapsed_sec, 3),
            "n_sites_in": len(result.sites_in) if result.sites_in else 0,
            "n_sites_out": len(result.sites_out) if result.sites_out else 0,
            "n_plots": len(result.plots),
            "outdir": str(result.outdir) if result.outdir else None,
            "steps": [
                {
                    "idx": sr.step_idx,
                    "label": sr.step_name,
                    "code": sr.step_code,
                    "ok": sr.ok,
                    "elapsed_sec": round(sr.elapsed_sec, 3),
                    "n_sites_in": sr.n_sites_in,
                    "n_sites_out": sr.n_sites_out,
                    "n_plots": len(sr.plots),
                    "error": str(sr.error) if sr.error else None,
                }
                for sr in result.step_results
            ],
        }
        return json.dumps(data, indent=2, default=str)

    if output_format == "csv":
        lines = [
            "idx,label,code,ok,elapsed_sec,n_sites_in,n_sites_out,n_plots,error"
        ]
        for sr in result.step_results:
            err = str(sr.error).replace(",", ";") if sr.error else ""
            lines.append(
                f"{sr.step_idx},{sr.step_name},{sr.step_code},"
                f"{sr.ok},{sr.elapsed_sec:.3f},"
                f"{sr.n_sites_in},{sr.n_sites_out},{len(sr.plots)},{err}"
            )
        return "\n".join(lines)

    # text
    return result.summary()


def _rich_pipe_table(
    title: str,
    rows: list[tuple[str, str]],
    style: str = "green",
) -> None:
    """Print a two-column key/value table (rich if available, plain fallback)."""
    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table import Table  # noqa: PLC0415

        console = Console()
        t = Table(title=title, border_style=style, show_header=False)
        t.add_column("Key", style="bold")
        t.add_column("Value", style="white")
        for k, v in rows:
            t.add_row(k, str(v))
        console.print(t)
    except ImportError:
        click.echo(f"\n{title}")
        for k, v in rows:
            click.echo(f"  {k:<30} {v}")


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


@click.group("pipe")
@click.pass_context
def pipe(ctx: click.Context) -> None:
    """Run, configure, and explore the pyCSAMT processing pipeline.

    \b
    Typical workflows:

      # Run a preset on the active survey
      pycsamt pipe run --preset full_processing

      # Run a config file on an explicit EDI directory
      pycsamt pipe run --config workflow.yaml --survey ./data/AMT/

      # Build and run a pipeline ad-hoc from the terminal
      pycsamt pipe run --steps NR001,FREQ002,FREQ001,FREQ004,SS001 --survey ./edis/

      # Validate a config without processing
      pycsamt pipe run --config workflow.yaml --dry-run

      # Explore what is available
      pycsamt pipe steps
      pycsamt pipe presets

      # Generate a starter config
      pycsamt pipe init --preset basic_qc --format yaml -o my_workflow.yaml
    """
    ctx.ensure_object(dict)
