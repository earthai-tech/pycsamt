# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli._base
=================

Root Click group.  All sub-commands are registered here and can reach
the live :data:`~pycsamt.api.cli.PYCSAMT_CLI` instance through
``ctx.obj["cli"]``.

The root group:

* reads ``~/.pycsamt.toml`` (if present) to set session defaults
* propagates ``--verbose`` / ``--no-color`` to :data:`PYCSAMT_CLI`
* prints a rich grouped help panel when invoked without a sub-command
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import click

from pycsamt import __version__

from ..api.cli.config import PYCSAMT_CLI, configure_cli


def _ensure_utf8_stdio() -> None:
    """Force UTF-8 on stdout/stderr.

    Windows falls back to the legacy console codepage (commonly cp1252)
    whenever stdout isn't attached to a real console — e.g. piped through
    ``| cat`` or redirected to a file — which crashes on the box-drawing
    characters used by the banner/help panel below.
    """
    for stream in (sys.stdout, sys.stderr):
        if stream is not None and hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except Exception:  # noqa: BLE001
                pass


_ensure_utf8_stdio()

# ---------------------------------------------------------------------------
# Optional rich formatting (degrades gracefully if rich is absent)
# ---------------------------------------------------------------------------
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text

    _RICH = True
except ImportError:
    _RICH = False


# ---------------------------------------------------------------------------
# Startup banner — block wordmark split at the "py" / "CSAMT" boundary and
# coloured with the site's brand palette (docs/source/_static/css/custom.css:
# --pycsamt-blue #3e65b0, --pst-color-primary #f15a29).
# ---------------------------------------------------------------------------
_BANNER_SPLIT = 17
_BANNER_LINES = (
    "██████╗ ██╗   ██╗ ██████╗███████╗ █████╗ ███╗   ███╗████████╗",
    "██╔══██╗╚██╗ ██╔╝██╔════╝██╔════╝██╔══██╗████╗ ████║╚══██╔══╝",
    "██████╔╝ ╚████╔╝ ██║     ███████╗███████║██╔████╔██║   ██║   ",
    "██╔═══╝   ╚██╔╝  ██║     ╚════██║██╔══██║██║╚██╔╝██║   ██║   ",
    "██║        ██║   ╚██████╗███████║██║  ██║██║ ╚═╝ ██║   ██║   ",
    "╚═╝        ╚═╝    ╚═════╝╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝   ╚═╝   ",
)
_BANNER_BLUE = "#3e65b0"
_BANNER_ORANGE = "#f15a29"


def _print_banner(console: "Console") -> None:
    """Print the pyCSAMT block wordmark, skipped on narrow terminals."""
    if console.size.width < len(_BANNER_LINES[0]):
        return
    for line in _BANNER_LINES:
        text = Text(line[:_BANNER_SPLIT], style=f"bold {_BANNER_BLUE}")
        text.append(line[_BANNER_SPLIT:].rstrip(), style=f"bold {_BANNER_ORANGE}")
        console.print(text)


def _print_rich_help(ctx: click.Context) -> None:
    """Print a grouped, coloured help panel using rich."""
    console = Console()
    _print_banner(console)

    table = Table.grid(padding=(0, 2))
    table.add_column(style="bold cyan", no_wrap=True)
    table.add_column()

    for name, cmd in ctx.command.commands.items():  # type: ignore[attr-defined]
        short_help = getattr(cmd, "help", "") or ""
        first_line = short_help.split("\n")[0].strip()
        table.add_row(name, first_line)

    panel = Panel(
        table,
        title=f"[bold]pyCSAMT {__version__}[/bold]",
        subtitle="[dim]pycsamt <command> --help  for command details[/dim]",
        border_style="bright_blue",
        expand=False,
    )
    console.print()
    console.print(panel)
    console.print()


def _load_toml_config(path: Path) -> dict[str, Any]:
    """Load ``~/.pycsamt.toml`` if it exists; return empty dict otherwise."""
    if not path.exists():
        return {}
    try:
        import tomllib  # Python 3.11+
    except ImportError:
        try:
            import tomli as tomllib  # type: ignore[no-redef]
        except ImportError:
            return {}
    with open(path, "rb") as fh:
        return tomllib.load(fh)


# ---------------------------------------------------------------------------
# Root group
# ---------------------------------------------------------------------------


class _RichGroup(click.Group):
    """Click Group that shows a rich help panel when called without args."""

    def invoke(self, ctx: click.Context) -> Any:
        if (
            not ctx.protected_args
            and not ctx.args
            and not ctx.invoked_subcommand
        ):
            if _RICH:
                _print_rich_help(ctx)
            else:
                click.echo(ctx.get_help())
            ctx.exit(0)
        return super().invoke(ctx)


@click.group(cls=_RichGroup, invoke_without_command=True)
@click.version_option(__version__, "--version", "-V", prog_name="pycsamt")
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity (-v info, -vv debug).",
)
@click.option(
    "--no-color",
    is_flag=True,
    default=False,
    help="Disable ANSI colour output.",
)
@click.pass_context
def main(ctx: click.Context, verbose: int, no_color: bool) -> None:
    """pyCSAMT — Audio-Magnetotelluric processing and inversion toolkit.

    Run  pycsamt <command> --help  for details on a specific command.
    """
    ctx.ensure_object(dict)

    # --- load ~/.pycsamt.toml and apply ALL sections (cli last so flags win) ---
    toml_path = Path.home() / ".pycsamt.toml"
    toml_cfg = _load_toml_config(toml_path)
    if toml_cfg:
        try:
            from pycsamt.cli.commands.config import (
                load_all_config,  # noqa: PLC0415
            )

            load_all_config(toml_cfg)
        except Exception:  # noqa: BLE001
            pass

    # CLI flags always win over config file
    if verbose:
        configure_cli(log__level=min(verbose, 2))
        from pycsamt.log.logger import (
            enable_console_logging,  # noqa: PLC0415
        )
        import logging  # noqa: PLC0415

        enable_console_logging(
            logging.DEBUG if verbose >= 2 else logging.INFO
        )
    if no_color:
        configure_cli(log__color=False)

    ctx.obj["cli"] = PYCSAMT_CLI


# ---------------------------------------------------------------------------
# Register sub-commands
# ---------------------------------------------------------------------------


def _register_commands() -> None:
    from pycsamt.cli.commands.avg import avg  # noqa: PLC0415
    from pycsamt.cli.commands.build import (
        build,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.config import (
        config,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.convert import (
        convert,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.edi import edi  # noqa: PLC0415
    from pycsamt.cli.commands.forward import (
        forward,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.info import (
        info,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.interp import (
        interp,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.invert import (
        invert,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.jones import (
        jones,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.map import (
        map,  # noqa: PLC0415,A001
    )
    from pycsamt.cli.commands.pipe import (
        pipe,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.site import (
        site,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.survey import (
        survey,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.tdem import (
        tdem,  # noqa: PLC0415
    )
    from pycsamt.cli.commands.transform import (
        transform,  # noqa: PLC0415
    )

    main.add_command(survey)
    main.add_command(site)
    main.add_command(pipe)
    main.add_command(map)
    main.add_command(info)
    main.add_command(convert)
    main.add_command(invert)
    main.add_command(transform)
    main.add_command(forward)
    main.add_command(interp)
    main.add_command(tdem)
    main.add_command(config)
    main.add_command(edi)
    main.add_command(jones)
    main.add_command(avg)
    main.add_command(build)


_register_commands()
