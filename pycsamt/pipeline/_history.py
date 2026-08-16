"""Append-only run history log for ``Pipeline.run(..., history=...)``.

One compact JSON line per run — not full ``StepResult`` objects — so the
log stays greppable and cheap to append/read even across many runs. This
is deliberately a plain JSON-lines file, not a database: it exists so a
user can answer "how did my last N runs compare," not to be a query engine.
"""

from __future__ import annotations

import datetime
import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ._pipeline import PipelineResult

__all__ = ["append_run", "load_history", "default_history_path"]


def default_history_path() -> Path:
    return Path.home() / ".pycsamt" / "pipeline_history.jsonl"


def _len_or_zero(obj: Any) -> int:
    try:
        return len(obj)
    except (TypeError, AttributeError):
        return 0


def _run_to_record(result: PipelineResult) -> dict:
    return {
        "timestamp": datetime.datetime.now(datetime.timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z"),
        "pipeline_name": result.pipeline_name,
        "ok": result.ok,
        "n_errors": result.n_errors,
        "elapsed_sec": round(result.elapsed_sec, 3),
        "n_sites_in": _len_or_zero(result.sites_in),
        "n_sites_out": _len_or_zero(result.sites_out),
        "outdir": str(result.outdir) if result.outdir is not None else None,
        "steps": [
            {
                "code": sr.step_code,
                "name": sr.step_name,
                "ok": sr.ok,
                "elapsed_sec": round(sr.elapsed_sec, 3),
                "cached": sr.cached,
            }
            for sr in result.step_results
        ],
    }


def append_run(path: str | Path, result: PipelineResult) -> None:
    """Append one JSON line summarizing *result* to *path*.

    Creates parent directories as needed. Never raises on a write failure
    when called from :meth:`Pipeline.run` — that call site wraps this in
    its own try/except and warns instead, the same "a logging feature must
    not abort the run" posture as every other optional side effect in the
    run loop (QC plots, intermediate snapshots).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    record = _run_to_record(result)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(record) + "\n")


def load_history(
    path: str | Path | None = None, *, last: int | None = None
) -> list[dict]:
    """Read back previously logged run summaries, oldest first.

    Parameters
    ----------
    path:
        Defaults to :func:`default_history_path`.
    last:
        When given, return only the most recent *last* entries.

    A malformed line (e.g. from a write interrupted mid-flush) is skipped
    rather than raising — the same "one bad entry is a non-event, not a
    crash" posture :class:`~pycsamt.pipeline.StepCache` already uses.
    """
    resolved = Path(path) if path is not None else default_history_path()
    if not resolved.exists():
        return []

    records: list[dict] = []
    with resolved.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue

    if last is not None:
        records = records[-last:]
    return records
