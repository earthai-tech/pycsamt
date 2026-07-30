# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.cli.commands.build._base
=================================

Root Click group for the external-solver build sub-commands, plus the
shared "locate a script, run it, forward its exit code" dispatcher used
by every leaf command in :mod:`pycsamt.cli.commands.build.commands`.

The actual compilation logic lives in the POSIX shell scripts shipped
under :mod:`pycsamt.models` (``models/_solver_build/*.sh``) -- see
:doc:`/user_guide/models/compilation`. This module does not duplicate
that logic; it only gives it a short, packaged, `pycsamt build ...`
entry point so a user never has to know or type the on-disk path to
those scripts.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import click

#: Windows locations checked for a Git-for-Windows bash.exe when one is
#: not already on PATH. Order matters: prefer the user-facing bash.exe
#: over the internal usr/bin one.
_WINDOWS_BASH_CANDIDATES = (
    r"C:\Program Files\Git\bin\bash.exe",
    r"C:\Program Files\Git\usr\bin\bash.exe",
    r"C:\Program Files (x86)\Git\bin\bash.exe",
)


def solver_build_dir() -> Path:
    """Return the on-disk path to ``pycsamt/models/_solver_build``.

    Resolved relative to the installed :mod:`pycsamt.models` package so
    this works identically from an editable checkout and from a built
    wheel (the directory is shipped via ``[tool.setuptools.package-data]``).
    """
    import pycsamt.models as _models  # noqa: PLC0415

    return Path(_models.__file__).resolve().parent / "_solver_build"


def find_bash() -> str | None:
    """Locate a bash-compatible shell to run the ``.sh`` build scripts."""
    for candidate in ("bash", "sh"):
        found = shutil.which(candidate)
        if found:
            return found
    if sys.platform == "win32":
        for guess in _WINDOWS_BASH_CANDIDATES:
            if Path(guess).exists():
                return guess
    return None


def run_solver_script(solver: str, extra_args: tuple[str, ...]) -> int:
    """Run ``_solver_build/<solver>.sh`` with *extra_args*, streaming output.

    Returns the child process's exit code so the caller can propagate it
    via ``sys.exit`` / ``raise SystemExit``.
    """
    script = solver_build_dir() / f"{solver}.sh"
    if not script.exists():
        click.echo(
            f"error: build script not found for '{solver}': {script}\n"
            "This pycsamt installation appears to be missing the "
            "models/_solver_build resources. Editable/dev checkouts and "
            "wheels built from the current pyproject.toml both include "
            "them; reinstalling pycsamt should fix this.",
            err=True,
        )
        return 1

    bash = find_bash()
    if bash is None:
        click.echo(
            "error: no bash-compatible shell found on PATH.\n"
            "  - Linux / macOS: bash is normally preinstalled.\n"
            "  - Windows: install Git for Windows "
            "(https://git-scm.com/download/win) and re-run this command "
            "from a shell that has it on PATH (e.g. Git Bash), or run it "
            "inside WSL.",
            err=True,
        )
        return 1

    cmd = [bash, str(script), *extra_args]
    try:
        return subprocess.call(cmd)
    except KeyboardInterrupt:
        return 130


@click.group("build")
@click.pass_context
def build(ctx: click.Context) -> None:
    """Compile the external EM solvers pyCSAMT integrates with.

    Thin, packaged front-end for the build scripts shipped under
    ``pycsamt/models/_solver_build/``. Each sub-command locates and runs
    the matching shell script for you, so you never need to know or type
    that on-disk path yourself. Every option after the sub-command name
    is forwarded verbatim to the underlying script -- run
    ``pycsamt build <solver> --help`` for the full list.

    \b
    Sub-commands:
      modem2d     Compile ModEM's 2-D MT inversion binary (Mod2DMT)
      modem3d     Compile ModEM's 3-D MT inversion binary (Mod3DMT)
      occam2d     Compile the Occam2D inversion binary
      mare2dem    Download / build MARE2DEM (Linux, macOS, or WSL only)

    \b
    Examples:
      pycsamt build modem3d --auto-install -y
      pycsamt build occam2d --clean
      pycsamt build mare2dem            # status only
      pycsamt build modem2d --help      # full option list for modem2d.sh
    """
    ctx.ensure_object(dict)
