# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.api.cli.params
======================

Custom Click parameter types for pyCSAMT-specific inputs.

Use these as the ``type=`` argument on ``@click.argument`` or
``@click.option`` to get automatic validation and helpful error messages.

Examples
--------
::

    import click
    from pycsamt.api.cli.params import EDIPath, FreqRange, StationList

    @click.command()
    @click.argument("edi", type=EDIPath())
    @click.option("--freq", type=FreqRange(), default=None)
    @click.option("--stations", type=StationList(), default=None)
    def my_command(edi, freq, stations):
        ...
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click

# ---------------------------------------------------------------------------
# EDIPath — single EDI file, must exist and have .edi extension
# ---------------------------------------------------------------------------


class EDIPath(click.Path):
    """A Click Path type that additionally requires a ``.edi`` extension."""

    name = "EDI_FILE"

    def __init__(self) -> None:
        super().__init__(
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            path_type=Path,
        )

    def convert(
        self,
        value: Any,
        param: click.Parameter | None,
        ctx: click.Context | None,
    ) -> Path:
        path: Path = super().convert(value, param, ctx)
        if path.suffix.lower() != ".edi":
            self.fail(
                f"{path.name!r} does not have a .edi extension.",
                param,
                ctx,
            )
        return path


# ---------------------------------------------------------------------------
# EDIDir — directory that contains at least one .edi file
# ---------------------------------------------------------------------------


class EDIDir(click.Path):
    """A Click Path type for a directory that contains EDI files."""

    name = "EDI_DIR"

    def __init__(self) -> None:
        super().__init__(
            exists=True,
            file_okay=False,
            dir_okay=True,
            readable=True,
            path_type=Path,
        )

    def convert(
        self,
        value: Any,
        param: click.Parameter | None,
        ctx: click.Context | None,
    ) -> Path:
        path: Path = super().convert(value, param, ctx)
        edi_files = list(path.glob("*.edi"))
        if not edi_files:
            self.fail(
                f"Directory {str(path)!r} contains no .edi files.",
                param,
                ctx,
            )
        return path


# ---------------------------------------------------------------------------
# FreqRange — "fmin:fmax" string → (float, float)
# ---------------------------------------------------------------------------


class FreqRange(click.ParamType):
    """Parse a ``fmin:fmax`` frequency-range string into ``(float, float)``.

    Examples
    --------
    ``--freq 0.1:1000``   →  ``(0.1, 1000.0)``
    """

    name = "FREQ_RANGE"

    def convert(
        self,
        value: Any,
        param: click.Parameter | None,
        ctx: click.Context | None,
    ) -> tuple[float, float]:
        if isinstance(value, tuple):
            return value
        try:
            lo_str, hi_str = str(value).split(":")
            lo, hi = float(lo_str), float(hi_str)
        except (ValueError, TypeError):
            self.fail(
                f"{value!r} is not a valid frequency range.  "
                "Expected format: fmin:fmax  (e.g. 0.1:1000)",
                param,
                ctx,
            )
        if lo >= hi:
            self.fail(
                f"fmin ({lo}) must be less than fmax ({hi}).",
                param,
                ctx,
            )
        return (lo, hi)


# ---------------------------------------------------------------------------
# StationList — comma-separated station ids → list[str]
# ---------------------------------------------------------------------------


class StationList(click.ParamType):
    """Parse a comma-separated list of station identifiers.

    Examples
    --------
    ``--stations S01,S02,S03``  →  ``['S01', 'S02', 'S03']``
    """

    name = "STATION_LIST"

    def convert(
        self,
        value: Any,
        param: click.Parameter | None,
        ctx: click.Context | None,
    ) -> list[str]:
        if isinstance(value, list):
            return value
        parts = [s.strip() for s in str(value).split(",") if s.strip()]
        if not parts:
            self.fail(
                f"{value!r} is not a valid station list.  "
                "Provide comma-separated station IDs, e.g. S01,S02,S03",
                param,
                ctx,
            )
        return parts


# ---------------------------------------------------------------------------
# PipeStepList — comma-separated step codes/names → list[str] of codes
# ---------------------------------------------------------------------------


class PipeStepList(click.ParamType):
    """Parse a comma-separated list of pipeline step codes or names.

    Each token is validated against the step registry, so the user gets
    an immediate error with a hint if they mistype a code.

    Examples
    --------
    ``--steps NR001,FREQ002,FREQ001,FREQ004,SS001``
        → ``['NR001', 'FREQ002', 'FREQ001', 'FREQ004', 'SS001']``

    ``--steps notch_powerline,drop_duplicates,align_grid``
        → ``['NR001', 'FREQ002', 'FREQ004']``
    """

    name = "STEP_LIST"

    def convert(
        self,
        value: Any,
        param: click.Parameter | None,
        ctx: click.Context | None,
    ) -> list[str]:
        if isinstance(value, list):
            return value
        tokens = [t.strip() for t in str(value).split(",") if t.strip()]
        if not tokens:
            self.fail(
                f"{value!r} is not a valid step list.  "
                "Provide comma-separated codes, e.g. NR001,FREQ001,SS001",
                param,
                ctx,
            )
        try:
            from pycsamt.pipeline import (
                lookup_step,  # noqa: PLC0415
            )
        except ImportError:
            return tokens  # best-effort if pipeline not yet importable
        codes: list[str] = []
        for token in tokens:
            try:
                spec = lookup_step(token)
                codes.append(spec.code)
            except KeyError:
                self.fail(
                    f"Unknown step {token!r}.  "
                    "Run  pycsamt pipe steps  to see all available steps.",
                    param,
                    ctx,
                )
        return codes


__all__ = [
    "EDIDir",
    "EDIPath",
    "FreqRange",
    "PipeStepList",
    "StationList",
]
