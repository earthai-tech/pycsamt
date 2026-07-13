# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""pycsamt invert status — inspect an inversion working directory."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import click

from ....api.cli.config import configure_cli
from ....api.cli.options import (
    format_option,
    no_color_option,
    verbose_option,
)
from ._base import _resolve_solver, _rich_table, invert


@invert.command("status")
@click.argument(
    "workdir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    metavar="WORKDIR",
)
@click.option(
    "--solver",
    type=click.Choice(["occam2d", "modem"], case_sensitive=False),
    default=None,
    help="Inversion code (auto-detected if omitted).",
)
@format_option
@verbose_option
@no_color_option
@click.pass_context
def status(
    ctx: click.Context,
    workdir: Path,
    solver: str | None,
    output_format: str,
    verbose: int,
    no_color: bool,
) -> None:
    """Inspect the state of an inversion working directory.

    Reports which input/output files are present and the current
    iteration count.  The solver is auto-detected from workdir contents.

    \b
    Examples:
      pycsamt invert status run01/
      pycsamt invert status run01/ --format json
      pycsamt invert status modem_run/ --solver modem
    """
    configure_cli(log__level=verbose, log__color=not no_color)
    solver_name = _resolve_solver(solver, workdir)
    info = _status_dict(workdir, solver_name, verbose)

    if output_format == "json":
        click.echo(json.dumps(info, indent=2, default=str))
    else:
        _print_status(info, solver_name)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _status_dict(workdir: Path, solver: str, verbose: int) -> dict[str, Any]:
    files = {p.name for p in workdir.iterdir() if p.is_file()}

    info: dict[str, Any] = {
        "workdir": str(workdir),
        "solver": solver,
        "files_present": sorted(files),
        "n_files": len(files),
    }

    if solver == "occam2d":
        iter_files = sorted(workdir.glob("*.iter"))
        resp_files = sorted(workdir.glob("*.resp"))
        info["n_iterations_done"] = len(iter_files)
        info["last_iter_file"] = iter_files[-1].name if iter_files else None
        info["n_resp_files"] = len(resp_files)
        info["has_data"] = bool(
            any(f.endswith(".dat") or "data" in f.lower() for f in files)
        )
        info["has_mesh"] = "Occam2DMesh" in files
        info["has_model"] = "Occam2DModel" in files
        info["has_startup"] = bool(any("startup" in f.lower() for f in files))
        info["ready_to_run"] = all(
            [
                info["has_data"],
                info["has_mesh"],
                info["has_model"],
                info["has_startup"],
            ]
        )
        log_path = workdir / "logFile.logFile"
        if not log_path.exists():
            logs = list(workdir.glob("*.log"))
            log_path = logs[-1] if logs else None
        info["rms_last"] = _last_rms_occam(log_path)

    else:
        model_files = sorted(workdir.glob("*.rho"))
        info["n_model_files"] = len(model_files)
        info["last_model_file"] = (
            model_files[-1].name if model_files else None
        )
        info["has_data"] = bool(any(".dat" in f for f in files))
        info["has_control"] = bool(any(".inv" in f for f in files))
        info["has_cov"] = bool(any(".cov" in f for f in files))
        info["ready_to_run"] = info["has_data"] and info["has_control"]
        info["rms_last"] = None

    return info


def _last_rms_occam(log_path: Path | None) -> float | None:
    if log_path is None or not log_path.exists():
        return None
    try:
        for line in reversed(
            log_path.read_text(errors="replace").splitlines()
        ):
            low = line.lower()
            if "rms" in low or "misfit" in low:
                for token in line.split():
                    try:
                        val = float(token)
                        if 0.0 < val < 1000.0:
                            return round(val, 4)
                    except ValueError:
                        pass
    except Exception:  # noqa: BLE001
        pass
    return None


def _print_status(info: dict[str, Any], solver: str) -> None:
    rows: list[tuple[str, str]] = [
        ("Workdir", info["workdir"]),
        ("Solver", solver.upper()),
        ("Files present", str(info["n_files"])),
    ]
    if solver == "occam2d":
        rows += [
            ("Has data file", "✓" if info["has_data"] else "✗"),
            ("Has mesh file", "✓" if info["has_mesh"] else "✗"),
            ("Has model file", "✓" if info["has_model"] else "✗"),
            ("Has startup file", "✓" if info["has_startup"] else "✗"),
            ("Ready to run", "✓" if info["ready_to_run"] else "✗"),
            ("Iterations done", str(info["n_iterations_done"])),
            ("Last .iter file", info["last_iter_file"] or "—"),
            ("Last RMS", str(info["rms_last"]) if info["rms_last"] else "—"),
        ]
    else:
        rows += [
            ("Has data file", "✓" if info["has_data"] else "✗"),
            ("Has control file", "✓" if info["has_control"] else "✗"),
            ("Has cov file", "✓" if info["has_cov"] else "✗"),
            ("Ready to run", "✓" if info["ready_to_run"] else "✗"),
            ("Model files", str(info["n_model_files"])),
            ("Last model file", info["last_model_file"] or "—"),
        ]
    _rich_table("Inversion Status", rows, style="blue")
