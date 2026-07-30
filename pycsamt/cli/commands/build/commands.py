# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.build.commands
====================================

Leaf commands for ``pycsamt build <solver>``.

Each one is a thin pass-through to
``pycsamt/models/_solver_build/<solver>.sh`` (see
:func:`pycsamt.cli.commands.build._base.run_solver_script`). Click's own
``--help`` handling is disabled on these commands on purpose: the shell
scripts already print a complete, hand-written help text (flags,
defaults, examples) via their own ``-h``/``--help``, and forwarding it
verbatim avoids maintaining the same option list in two places.
"""

from __future__ import annotations

import click

from ._base import build, run_solver_script

# Unknown options (e.g. --auto-install, --mpi) must reach the script
# unparsed, and --help/-h must NOT be swallowed by click -- both are
# forwarded to the shell script itself.
_CONTEXT_SETTINGS = {"ignore_unknown_options": True}

_SOLVERS = {
    "modem2d": "Compile ModEM's 2-D MT inversion binary (Mod2DMT).",
    "modem3d": "Compile ModEM's 3-D MT inversion binary (Mod3DMT).",
    "occam2d": "Compile the Occam2D inversion binary.",
    "mare2dem": "Download / build MARE2DEM (Linux, macOS, or WSL only).",
}


def _make_command(name: str, summary: str) -> click.Command:
    @click.command(
        name,
        context_settings=_CONTEXT_SETTINGS,
        add_help_option=False,
        help=(
            f"{summary}\n\n"
            f"All options are forwarded to _solver_build/{name}.sh -- run "
            f"'pycsamt build {name} --help' for the full option list, "
            "including --auto-install, --clean, and platform-specific "
            "flags."
        ),
    )
    @click.argument("args", nargs=-1, type=click.UNPROCESSED)
    def _cmd(args: tuple[str, ...]) -> None:
        raise SystemExit(run_solver_script(name, args))

    return _cmd


modem2d = _make_command("modem2d", _SOLVERS["modem2d"])
modem3d = _make_command("modem3d", _SOLVERS["modem3d"])
occam2d = _make_command("occam2d", _SOLVERS["occam2d"])
mare2dem = _make_command("mare2dem", _SOLVERS["mare2dem"])

for _cmd in (modem2d, modem3d, occam2d, mare2dem):
    build.add_command(_cmd)
