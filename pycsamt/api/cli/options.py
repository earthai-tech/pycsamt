# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.api.cli.options
=======================

Reusable Click option decorators shared across all pyCSAMT commands.

Importing a decorator and stacking it onto a Click command automatically
wires up the matching argument name — no need to repeat option strings in
every command.

Usage
-----
::

    import click
    from pycsamt.api.cli.options import verbose_option, output_dir_option


    @click.command()
    @verbose_option
    @output_dir_option
    def my_command(verbose, output_dir): ...
"""

from __future__ import annotations

from pathlib import Path

import click

# ---------------------------------------------------------------------------
# Logging / display options
# ---------------------------------------------------------------------------

verbose_option = click.option(
    "-v",
    "--verbose",
    count=True,
    help=(
        "Increase output verbosity.  "
        "Use -v for informational messages, -vv for debug output."
    ),
)

no_color_option = click.option(
    "--no-color",
    is_flag=True,
    default=False,
    help="Disable ANSI colour in terminal output.",
)

# ---------------------------------------------------------------------------
# Output options
# ---------------------------------------------------------------------------

output_dir_option = click.option(
    "-o",
    "--output-dir",
    type=click.Path(file_okay=False, writable=True, path_type=Path),
    default=Path("."),
    show_default=True,
    help="Directory where output files are written.",
)

format_option = click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["text", "json", "csv"], case_sensitive=False),
    default="text",
    show_default=True,
    help="Output format for command results.",
)

overwrite_option = click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    help="Overwrite existing output files without prompting.",
)

# ---------------------------------------------------------------------------
# Survey / data-source options
# ---------------------------------------------------------------------------

survey_option = click.option(
    "-S",
    "--survey",
    "survey_path",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    metavar="EDI_DIR",
    help=(
        "EDI directory or file to use as the survey source.  "
        "When omitted the active survey set via "
        "'pycsamt survey set' is used automatically."
    ),
)

fresh_option = click.option(
    "--fresh",
    is_flag=True,
    default=False,
    help=(
        "Force re-parse of the survey from disk, ignoring the cache.  "
        "Useful after editing EDI files."
    ),
)

# ---------------------------------------------------------------------------
# Processing / build options
# ---------------------------------------------------------------------------

n_jobs_option = click.option(
    "-j",
    "--jobs",
    "n_jobs",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    help="Number of parallel worker processes.",
)

no_cache_option = click.option(
    "--no-cache",
    is_flag=True,
    default=False,
    help="Disable intermediate result caching.",
)

# ``click.option`` returns decorator callables without useful docstrings.
# Assign concise API descriptions so autosummary can document these public
# convenience decorators like ordinary functions.
verbose_option.__doc__ = (
    "Add the repeatable ``--verbose`` logging option to a Click command."
)
no_color_option.__doc__ = (
    "Add the ``--no-color`` terminal-output option to a Click command."
)
output_dir_option.__doc__ = (
    "Add the writable ``--output-dir`` path option to a Click command."
)
format_option.__doc__ = (
    "Add the text, JSON, or CSV output-format option to a Click command."
)
overwrite_option.__doc__ = (
    "Add the ``--overwrite`` confirmation-bypass option to a Click command."
)
n_jobs_option.__doc__ = (
    "Add the positive ``--jobs`` parallel-worker option to a Click command."
)
no_cache_option.__doc__ = "Add the ``--no-cache`` processing option to a Click command."

# ---------------------------------------------------------------------------
# Convenience bundle: options attached to nearly every command
# ---------------------------------------------------------------------------


def common_options(f):
    """Attach verbose, no-color, format, and output-dir to a command."""
    for decorator in reversed(
        [
            verbose_option,
            no_color_option,
            format_option,
            output_dir_option,
        ]
    ):
        f = decorator(f)
    return f


__all__ = [
    "common_options",
    "format_option",
    "fresh_option",
    "n_jobs_option",
    "no_cache_option",
    "no_color_option",
    "output_dir_option",
    "overwrite_option",
    "survey_option",
    "verbose_option",
]
