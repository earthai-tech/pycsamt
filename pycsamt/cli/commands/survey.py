# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt survey
==============

Manage the active survey context — the EDI collection / Sites object
that all other commands use by default when no explicit path is given.

Sub-commands
------------
set       Register an EDI directory (or file) as the active survey and
          pre-build the Sites cache.
show      Print a summary of the currently active survey.
clear     Unset the active survey (cache files are preserved).
rebuild   Force-rebuild the Sites cache for the current active survey.
cache     Low-level cache management (list, purge).

Usage
-----
::

    # Point the CLI at a survey once
    pycsamt survey set data/willy/

    # Every subsequent command picks it up automatically
    pycsamt info
    pycsamt invert build --solver occam2d --workdir run01/

    # Override for a single command
    pycsamt invert build --survey data/other_survey/ --workdir run02/

    # Inspect active context
    pycsamt survey show
    pycsamt survey show --format json

    # Invalidate / rebuild
    pycsamt survey rebuild
    pycsamt survey rebuild --force    # unconditional

    # Remove active context (cache is kept)
    pycsamt survey clear

    # Cache management
    pycsamt survey cache list
    pycsamt survey cache purge        # purge active survey cache only
    pycsamt survey cache purge --all  # purge every cached survey
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import click

from ...api.cli.config import configure_cli
from ...api.cli.options import (
    format_option,
    no_color_option,
    verbose_option,
)
from ..survey import (
    _CACHE_ROOT,
    SurveyContext,
    _purge_cache,
    set_survey,
    survey_summary,
)

# ---------------------------------------------------------------------------
# Rich helpers
# ---------------------------------------------------------------------------


def _rich_kv(
    title: str, rows: list[tuple[str, str]], style: str = "cyan"
) -> None:
    try:
        from rich.console import Console  # noqa: PLC0415
        from rich.table import Table  # noqa: PLC0415

        console = Console()
        t = Table(
            title=title,
            border_style=style,
            show_header=False,
            box=None,
            padding=(0, 2),
        )
        t.add_column(style="bold dim")
        t.add_column(style="white")
        for k, v in rows:
            t.add_row(k, v)
        console.print()
        console.print(t)
        console.print()
    except ImportError:
        click.echo(f"\n{title}")
        for k, v in rows:
            click.echo(f"  {k:<26} {v}")
        click.echo()


def _rich_warn(msg: str) -> None:
    try:
        from rich.console import Console  # noqa: PLC0415

        Console().print(f"[yellow]Warning:[/yellow] {msg}")
    except ImportError:
        click.echo(f"Warning: {msg}", err=True)


def _rich_ok(msg: str) -> None:
    try:
        from rich.console import Console  # noqa: PLC0415

        Console().print(f"[green]✓[/green] {msg}")
    except ImportError:
        click.echo(msg)


# ---------------------------------------------------------------------------
# survey — top-level group
# ---------------------------------------------------------------------------


@click.group("survey")
@click.pass_context
def survey(ctx: click.Context) -> None:
    """Manage the active survey (EDI collection) used by all commands.

    \b
    Quick start:
      pycsamt survey set  data/willy/     # register + cache
      pycsamt survey show                 # inspect
      pycsamt survey rebuild              # refresh cache
      pycsamt survey clear                # unset
    """
    ctx.ensure_object(dict)


# ---------------------------------------------------------------------------
# survey set
# ---------------------------------------------------------------------------


@survey.command("set")
@click.argument(
    "path", type=click.Path(exists=True, path_type=Path), metavar="EDI_DIR"
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Rebuild cache even if it is still valid.",
)
@verbose_option
@no_color_option
@click.pass_context
def survey_set(
    ctx: click.Context,
    path: Path,
    force: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Register EDI_DIR as the active survey and pre-build the Sites cache.

    After this command any pycsamt command that requires EDI data will
    use this survey automatically unless you pass an explicit path.

    \b
    Examples:
      pycsamt survey set data/willy/
      pycsamt survey set data/willy/ --force     # unconditional rebuild
      pycsamt survey set data/willy/ -v          # verbose build progress
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    click.echo(f"Setting active survey: {path.resolve()}")

    sites = set_survey(path, force=force, verbose=verbose)

    try:
        station_list = ", ".join(s.name for s in sites)
    except Exception:  # noqa: BLE001
        station_list = "—"

    _rich_ok(f"Active survey set — {len(sites)} station(s): {station_list}")


# ---------------------------------------------------------------------------
# survey show
# ---------------------------------------------------------------------------


@survey.command("show")
@format_option
@no_color_option
@click.pass_context
def survey_show(
    ctx: click.Context,
    output_format: str,
    no_color: bool,
) -> None:
    """Print a summary of the currently active survey.

    \b
    Examples:
      pycsamt survey show
      pycsamt survey show --format json
    """
    configure_cli(log__color=not no_color)
    summary = survey_summary()

    if summary is None:
        click.echo(
            "No active survey is set.\n\n  Run:  pycsamt survey set <edi_dir>",
            err=False,
        )
        return

    if output_format == "json":
        click.echo(json.dumps(summary, indent=2, default=str))
        return

    cache_status = (
        "[green]valid[/green]"
        if summary["cache_valid"]
        else "[red]stale[/red]"
    )
    try:
        from rich.markup import escape  # noqa: PLC0415

        escape(summary["survey_path"])
    except ImportError:
        summary["survey_path"]
        cache_status = "valid" if summary["cache_valid"] else "stale"

    rows = [
        ("Path", summary["survey_path"]),
        ("Set at", summary["set_at"]),
        ("Stations", str(summary["n_stations"])),
        ("Station list", ", ".join(summary["station_names"]) or "—"),
        ("Cache", cache_status),
        ("Cache path", summary["cache_path"]),
    ]
    if "cached_at" in summary:
        rows.append(("Cached at", summary["cached_at"]))

    _rich_kv("Active Survey", rows, style="bright_cyan")


# ---------------------------------------------------------------------------
# survey clear
# ---------------------------------------------------------------------------


@survey.command("clear")
@click.option(
    "--yes", is_flag=True, default=False, help="Skip the confirmation prompt."
)
@no_color_option
@click.pass_context
def survey_clear(
    ctx: click.Context,
    yes: bool,
    no_color: bool,
) -> None:
    """Unset the active survey context.

    Cache files are preserved.  Run ``pycsamt survey cache purge`` to
    also delete the cached Sites objects.

    \b
    Examples:
      pycsamt survey clear
      pycsamt survey clear --yes   # no prompt
    """
    configure_cli(log__color=not no_color)
    active = SurveyContext.load()
    if active is None:
        click.echo("No active survey is set — nothing to clear.")
        return

    if not yes:
        confirmed = click.confirm(
            f"Clear active survey {active.survey_path}?", default=False
        )
        if not confirmed:
            click.echo("Aborted.")
            return

    SurveyContext.clear()
    _rich_ok("Active survey cleared.")


# ---------------------------------------------------------------------------
# survey rebuild
# ---------------------------------------------------------------------------


@survey.command("rebuild")
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Rebuild unconditionally even if cache is still valid.",
)
@verbose_option
@no_color_option
@click.pass_context
def survey_rebuild(
    ctx: click.Context,
    force: bool,
    verbose: int,
    no_color: bool,
) -> None:
    """Rebuild the Sites cache for the active survey.

    Use this after editing EDI files, adding stations, or when
    the cache appears stale.

    \b
    Examples:
      pycsamt survey rebuild
      pycsamt survey rebuild --force   # unconditional
      pycsamt survey rebuild -v        # show progress
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    active = SurveyContext.load()
    if active is None:
        raise click.UsageError(
            "No active survey to rebuild.\n"
            "  Run:  pycsamt survey set <edi_dir>"
        )

    click.echo(f"Rebuilding cache for {active.survey_path} …")
    sites = set_survey(active.survey_path, force=True, verbose=verbose)
    _rich_ok(f"Cache rebuilt — {len(sites)} station(s).")


# ---------------------------------------------------------------------------
# survey cache (sub-group)
# ---------------------------------------------------------------------------


@survey.group("cache")
@click.pass_context
def survey_cache(ctx: click.Context) -> None:
    """Low-level cache management.

    \b
    Sub-commands:
      list    Show all cached surveys.
      purge   Delete cache for the active survey (or all surveys).
    """
    ctx.ensure_object(dict)


@survey_cache.command("list")
@format_option
@click.pass_context
def cache_list(ctx: click.Context, output_format: str) -> None:
    """List all cached survey entries.

    \b
    Examples:
      pycsamt survey cache list
      pycsamt survey cache list --format json
    """
    if not _CACHE_ROOT.exists():
        click.echo("No cache directory found.")
        return

    entries = []
    for d in sorted(_CACHE_ROOT.iterdir()):
        if not d.is_dir():
            continue
        meta_path = d / "meta.json"
        pkl_path = d / "sites.pkl"
        entry: dict = {"key": d.name}
        if meta_path.exists():
            try:
                m = json.loads(meta_path.read_text(encoding="utf-8"))
                entry.update(m)
            except Exception:  # noqa: BLE001
                pass
        entry["pkl_size_kb"] = (
            round(pkl_path.stat().st_size / 1024, 1)
            if pkl_path.exists()
            else 0
        )
        entries.append(entry)

    if not entries:
        click.echo("Cache is empty.")
        return

    if output_format == "json":
        click.echo(json.dumps(entries, indent=2, default=str))
        return

    # text table
    active_ctx = SurveyContext.load()
    active_key = active_ctx.cache_key if active_ctx else None

    rows = []
    for e in entries:
        marker = " [active]" if e["key"] == active_key else ""
        path = e.get("survey_path", "—")
        n = str(e.get("n_stations", "?"))
        size = f"{e['pkl_size_kb']} KB"
        cached = e.get("cached_at", "—")
        rows.append(
            (e["key"] + marker, f"{path}  ({n} stations, {size}, {cached})")
        )

    _rich_kv("Cached Surveys", rows, style="blue")


@survey_cache.command("purge")
@click.option(
    "--all",
    "purge_all",
    is_flag=True,
    default=False,
    help="Purge every cached survey, not just the active one.",
)
@click.option(
    "--yes", is_flag=True, default=False, help="Skip the confirmation prompt."
)
@click.pass_context
def cache_purge(ctx: click.Context, purge_all: bool, yes: bool) -> None:
    """Delete the Sites cache for the active survey (or all surveys).

    \b
    Examples:
      pycsamt survey cache purge          # active survey only
      pycsamt survey cache purge --all    # every cached survey
      pycsamt survey cache purge --yes    # no confirmation prompt
    """
    if purge_all:
        if not _CACHE_ROOT.exists():
            click.echo("Cache directory does not exist — nothing to purge.")
            return
        if not yes:
            confirmed = click.confirm(
                f"Delete ALL cached surveys in {_CACHE_ROOT}?", default=False
            )
            if not confirmed:
                click.echo("Aborted.")
                return
        shutil.rmtree(_CACHE_ROOT)
        _rich_ok("All cached surveys purged.")
        return

    active = SurveyContext.load()
    if active is None:
        raise click.UsageError(
            "No active survey.  Pass --all to purge every cache entry."
        )
    if not yes:
        confirmed = click.confirm(
            f"Delete cache for {active.survey_path}?", default=False
        )
        if not confirmed:
            click.echo("Aborted.")
            return
    _purge_cache(active.cache_key)
    _rich_ok(f"Cache purged for {active.survey_path}.")
