# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt format validate
=======================

Structural + round-trip check of a ``.pcsf`` / ``.pcsm`` file:

1. header peek (``peek_kind``);
2. full load (``read_pcsf_or_pcsm``);
3. schema check (``PCSFModel.validate``);
4. round-trip: re-write to a temp file of the same encoding, re-read,
   and compare the resistivity array element-for-element.

Exits non-zero if any step fails.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import click

from ._base import _rich_table, fmt


@fmt.command("validate")
@click.argument(
    "file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    metavar="FILE",
)
@click.option(
    "--no-roundtrip",
    is_flag=True,
    default=False,
    help="Skip the write → read → compare round-trip step.",
)
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["text", "json"], case_sensitive=False),
    default="text",
    show_default=True,
)
def validate(file: Path, no_roundtrip: bool, output_format: str) -> None:
    """Validate a PCSF / PCSM file.

    \b
    Examples:
      pycsamt format validate model.pcsf
      pycsamt format validate model.pcsm --no-roundtrip
    """
    import numpy as np

    from pycsamt.format.io import write_pcsf
    from pycsamt.format.text import (
        peek_kind,
        read_pcsf_or_pcsm,
        write_pcsm,
    )

    name = file.name.lower()
    is_pcsm = name.endswith((".pcsm", ".pcsm.gz"))
    if not (is_pcsm or name.endswith(".pcsf")):
        raise click.ClickException(
            f"{file.name} is not a .pcsf / .pcsm / .pcsm.gz file."
        )

    checks: list[dict] = []

    def _record(step: str, ok: bool, detail: str = "") -> bool:
        checks.append({"step": step, "ok": ok, "detail": detail})
        return ok

    kind = None
    try:
        kind = peek_kind(file)
        _record("header", True, f"geometry={kind}")
    except Exception as exc:  # noqa: BLE001
        _record("header", False, str(exc))

    model = None
    if kind is not None:
        try:
            model = read_pcsf_or_pcsm(file)
            _record("load", True, f"kind={model.kind}")
        except Exception as exc:  # noqa: BLE001
            _record("load", False, str(exc))

    if model is not None:
        try:
            model.validate()
            _record("schema", True)
        except Exception as exc:  # noqa: BLE001
            _record("schema", False, str(exc))

    if model is not None and not no_roundtrip:
        try:
            with tempfile.TemporaryDirectory() as tmp:
                if is_pcsm:
                    rt = Path(tmp) / "rt.pcsm"
                    write_pcsm(model, rt)
                else:
                    rt = Path(tmp) / "rt.pcsf"
                    write_pcsf(model, rt)
                back = read_pcsf_or_pcsm(rt)
                a, b = model.resistivity, back.resistivity
                if a is None and b is None:
                    _record("roundtrip", True, "no resistivity array")
                elif a is None or b is None:
                    _record("roundtrip", False, "resistivity presence changed")
                else:
                    np.testing.assert_allclose(
                        np.asarray(a), np.asarray(b), rtol=1e-9, atol=0,
                        equal_nan=True,
                    )
                    _record(
                        "roundtrip",
                        True,
                        f"resistivity {list(a.shape)} matches",
                    )
        except Exception as exc:  # noqa: BLE001
            _record("roundtrip", False, str(exc))

    ok = all(c["ok"] for c in checks) and len(checks) > 0
    result = {"file": str(file), "valid": ok, "checks": checks}

    if output_format == "json":
        click.echo(json.dumps(result, indent=2, default=str))
    else:
        _rich_table(
            "format validate",
            [
                (
                    ("✓ " if c["ok"] else "✗ ") + c["step"],
                    c["detail"] or ("pass" if c["ok"] else "FAIL"),
                )
                for c in checks
            ],
            style="green" if ok else "red",
        )
        click.echo(
            f"\n{'✓ VALID' if ok else '✗ INVALID'} — {file}"
        )

    if not ok:
        raise SystemExit(1)
