# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt pipe run — run a processing pipeline against MT/AMT EDI data."""

from __future__ import annotations

import sys
from pathlib import Path

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    fresh_option,
    n_jobs_option,
    no_color_option,
    survey_option,
    verbose_option,
)
from ....api.cli.params import PipeStepList
from ._base import (
    _format_run_result,
    _resolve_pipeline,
    _resolve_sites,
    _trim_pipeline,
    pipe,
)


@pipe.command("run")
# ── EDI source (positional + flag) ────────────────────────────────────────
@click.argument(
    "edi_source",
    type=click.Path(exists=True, path_type=Path),
    metavar="[EDI_DIR]",
    required=False,
    default=None,
)
@survey_option
@fresh_option
# ── Pipeline definition ───────────────────────────────────────────────────
@click.option(
    "--config",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    metavar="FILE",
    help=(
        "Pipeline config file (.yaml, .json, or .py).  "
        "Takes precedence over --preset and --steps."
    ),
)
@click.option(
    "--preset",
    default=None,
    metavar="NAME",
    help=(
        "Named preset (basic_qc, noise_reduction, full_processing, "
        "tensor_analysis, dimensionality_filter, publication_ready, "
        "stratagem_mt, mt_qc, amt_qc, csamt_qc, csumt_qc).  "
        "Run  pycsamt pipe presets  to list all."
    ),
)
@click.option(
    "--steps",
    type=PipeStepList(),
    default=None,
    metavar="CODE,CODE,…",
    help=(
        "Comma-separated step codes to execute ad-hoc.  "
        "Example: NR001,FREQ002,FREQ001,FREQ004,SS001.  "
        "Run  pycsamt pipe steps  to list all codes."
    ),
)
@click.option(
    "--name",
    default=None,
    metavar="TEXT",
    help="Override the pipeline name stored in the output YAML.",
)
# ── Execution controls ────────────────────────────────────────────────────
@click.option(
    "--from-step",
    default=None,
    metavar="LABEL_OR_CODE",
    help="Start execution from this step (skip all earlier steps).",
)
@click.option(
    "--until-step",
    default=None,
    metavar="LABEL_OR_CODE",
    help="Stop execution after this step (discard all later steps).",
)
@click.option(
    "--n-steps",
    type=click.IntRange(min=1),
    default=None,
    metavar="INT",
    help="Run only the first N steps (useful for incremental debugging).",
)
@click.option(
    "--dry-run",
    is_flag=True,
    default=False,
    help=(
        "Resolve and display the pipeline and site count without "
        "running any processing.  No files are written."
    ),
)
@click.option(
    "--on-error",
    type=click.Choice(["raise", "warn", "skip"], case_sensitive=False),
    default="warn",
    show_default=True,
    help="Behaviour when a step raises an error.",
)
# ── Output controls ───────────────────────────────────────────────────────
@click.option(
    "--out",
    "out_dir",
    type=click.Path(path_type=Path),
    default=None,
    metavar="DIR",
    help=(
        "Root output directory for processed EDIs, plots, and reports.  "
        "Defaults to the pipeline config value or pipe_results/."
    ),
)
@click.option(
    "--no-plots",
    is_flag=True,
    default=False,
    help="Skip generation and saving of QC figures.",
)
@click.option(
    "--no-edi",
    is_flag=True,
    default=False,
    help="Skip writing processed EDI files to <out>/processed/.",
)
@click.option(
    "--no-report",
    is_flag=True,
    default=False,
    help="Skip writing HTML and text run reports.",
)
@click.option(
    "--dpi",
    type=click.IntRange(min=50),
    default=150,
    show_default=True,
    metavar="INT",
    help="DPI for saved QC figures.",
)
@click.option(
    "--plot-fmt",
    type=click.Choice(["png", "pdf", "svg"], case_sensitive=False),
    default="png",
    show_default=True,
    help="File format for saved QC figures.",
)
@click.option(
    "--cache",
    is_flag=True,
    default=False,
    help=(
        "Cache each step's output, keyed by the exact upstream data + step "
        "code + params.  A rerun of the identical command replays already-"
        "completed steps from cache instead of recomputing them — this is "
        "also how an interrupted run resumes.  Off by default: a step whose "
        "function is non-deterministic or reads state outside its own "
        "params is not safe to cache."
    ),
)
@click.option(
    "--cache-dir",
    type=click.Path(path_type=Path),
    default=None,
    metavar="DIR",
    help="Cache location.  Defaults to ~/.pycsamt/pipeline_cache.  Implies --cache.",
)
@click.option(
    "--live",
    is_flag=True,
    default=False,
    help=(
        "Render a live-updating status table (pending/running/OK/ERR/cached "
        "per step) instead of a static progress bar.  Implies visible "
        "progress even without -v."
    ),
)
@click.option(
    "--history",
    is_flag=True,
    default=False,
    help=(
        "Append a one-line JSON summary of this run to the run history log "
        "(default ~/.pycsamt/pipeline_history.jsonl).  See  pycsamt pipe history."
    ),
)
@click.option(
    "--history-file",
    type=click.Path(path_type=Path),
    default=None,
    metavar="FILE",
    help="Run-history log location.  Implies --history.",
)
@click.option(
    "--dashboard",
    is_flag=True,
    default=False,
    help=(
        "Also write dashboard.html — a richer, branded report with KPI "
        "stat tiles and inline-SVG charts (step status, per-step "
        "duration, site-count flow) — alongside the default report.html "
        "and summary.txt."
    ),
)
# ── Standard global options ───────────────────────────────────────────────
@n_jobs_option
@format_option
@verbose_option
@no_color_option
@click.pass_context
def run(
    ctx: click.Context,
    edi_source: Path | None,
    survey_path: Path | None,
    fresh: bool,
    config: Path | None,
    preset: str | None,
    steps: list[str] | None,
    name: str | None,
    from_step: str | None,
    until_step: str | None,
    n_steps: int | None,
    dry_run: bool,
    on_error: str,
    out_dir: Path | None,
    no_plots: bool,
    no_edi: bool,
    no_report: bool,
    dpi: int,
    plot_fmt: str,
    cache: bool,
    cache_dir: Path | None,
    live: bool,
    history: bool,
    history_file: Path | None,
    dashboard: bool,
    n_jobs: int,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Run a processing pipeline against MT/AMT sites.

    \b
    EDI source (highest → lowest priority):
      Positional EDI_DIR  →  --survey DIR  →  active survey context

    Pipeline definition (highest → lowest priority):
      --config FILE  →  --preset NAME  →  --steps CODE,CODE,...

    \b
    Examples:
      # Preset on the active survey (set with 'pycsamt survey set ./edis/')
      pycsamt pipe run --preset full_processing

      # Config file on an explicit directory
      pycsamt pipe run --config workflow.yaml --survey ./data/AMT/ --out results/

      # Fully ad-hoc pipeline from the terminal
      pycsamt pipe run --steps NR001,FREQ002,FREQ001,FREQ004,SS001 \\
              --survey ./edis/ --out results/

      # Validate without running anything
      pycsamt pipe run --config workflow.yaml --dry-run

      # Debug: run only the first 3 steps
      pycsamt pipe run --config workflow.yaml --n-steps 3 --dry-run

      # Partial run: skip early steps, stop at static-shift correction
      pycsamt pipe run --config workflow.yaml --from-step align --until-step correct_ss

      # Override error policy and output DPI
      pycsamt pipe run --preset basic_qc --on-error raise --dpi 300 --plot-fmt pdf

      # Cache step outputs — a rerun resumes instead of recomputing
      pycsamt pipe run --config workflow.yaml --survey ./edis/ --cache

      # Live per-step status table instead of a static progress bar
      pycsamt pipe run --preset basic_qc --live

      # Log this run, then compare later with  pycsamt pipe history
      pycsamt pipe run --preset basic_qc --history

      # Also write a branded dashboard.html with stat tiles and charts
      pycsamt pipe run --preset basic_qc --dashboard
    """
    configure_cli(log__level=verbose, log__color=not no_color)

    from pycsamt.api.pipe import (
        configure_pipe,  # noqa: PLC0415
    )

    # Apply pipeline-level config overrides
    pipe_overrides = dict(
        on_step_error=on_error,
        plot_dpi=dpi,
        plot_fmt=plot_fmt.lower(),
        show_progress=(verbose >= 1) or live,
    )
    if live:
        # --live is an explicit request to see it — don't disturb the
        # user's own progress_style configuration otherwise.
        pipe_overrides["progress_style"] = "rich"
    if dashboard:
        # --dashboard adds the "dashboard" report format for this run only
        # — don't disturb any other formats already configured.
        pipe_overrides["report_formats"] = ("html", "txt", "dashboard")
    configure_pipe(**pipe_overrides)

    # ── 1. Resolve pipeline ───────────────────────────────────────────────
    try:
        pipeline = _resolve_pipeline(config, preset, steps, name, verbose)
    except click.UsageError:
        raise
    except Exception as exc:
        raise click.ClickException(f"Failed to load pipeline: {exc}") from exc

    # ── 2. Apply step slicing ─────────────────────────────────────────────
    try:
        pipeline = _trim_pipeline(pipeline, from_step, until_step, n_steps)
    except click.BadParameter:
        raise

    if len(pipeline) == 0:
        raise click.ClickException(
            "The resolved pipeline has no steps after applying "
            "--from-step / --until-step / --n-steps."
        )

    # ── 3. Resolve sites ─────────────────────────────────────────────────
    try:
        sites = _resolve_sites(edi_source, survey_path, fresh, verbose)
    except click.UsageError:
        raise
    except Exception as exc:
        raise click.ClickException(f"Failed to load sites: {exc}") from exc

    n_sites = len(sites) if hasattr(sites, "__len__") else "?"
    cache_arg: bool | Path = cache_dir if cache_dir is not None else cache
    history_arg: bool | Path = history_file if history_file is not None else history

    # ── 4. Dry-run: print resolved pipeline and exit ──────────────────────
    if dry_run:
        click.echo(str(pipeline))
        click.echo()
        click.echo(f"Sites   : {n_sites}")
        click.echo(f"Steps   : {len(pipeline)}")
        click.echo(f"Out dir : {out_dir or 'pipe_results (default)'}")
        click.echo(f"Cache   : {cache_arg if cache_arg else 'disabled'}")
        click.echo(f"History : {history_arg if history_arg else 'disabled'}")
        click.echo(f"Dashboard : {'enabled' if dashboard else 'disabled'}")
        click.echo()
        click.echo("Dry run — no processing performed.")
        return

    # ── 5. Run ────────────────────────────────────────────────────────────
    if verbose >= 1:
        click.echo(f"  Processing {n_sites} site(s)…", err=True)

    try:
        result = pipeline.run(
            sites,
            outdir=out_dir,
            save_plots=not no_plots,
            save_edis=not no_edi,
            save_report=not no_report,
            cache=cache_arg,
            history=history_arg,
        )
    except Exception as exc:
        raise click.ClickException(f"Pipeline run failed: {exc}") from exc

    # ── 6. Output ─────────────────────────────────────────────────────────
    if not result.ok and verbose >= 1:
        click.echo(
            f"  Warning: {result.n_errors} step(s) raised errors.", err=True
        )

    summary_text = _format_run_result(result, output_format)
    click.echo(summary_text)

    if result.outdir and verbose >= 1:
        click.echo(f"\n  Output → {result.outdir}/", err=True)

    if not result.ok:
        sys.exit(1)
