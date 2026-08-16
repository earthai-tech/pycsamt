"""Live rich.table.Table rendering for ``progress_style="rich"``.

Kept separate from ``_pipeline.py`` so the run loop only has to call
:func:`make_rich_live` once and mutate the returned row dicts — the table
itself is always rebuilt from scratch on each repaint (``rich.table.Table``
has no supported in-place cell-update API), which is the standard idiom for
driving a changing table through :class:`rich.live.Live`.
"""

from __future__ import annotations

from typing import Any

__all__ = ["make_rich_live"]

_STATUS_STYLE = {
    "pending": None,
    "running": "yellow",
    "OK": "bold green",
    "ERR": "bold red",
}


def make_rich_live(pipeline_name: str, steps: list[tuple[str, Any]]):
    """Return ``(live, rows, render)`` for a pipeline about to run *steps*.

    *rows* is a list of plain dicts, one per step, pre-populated as
    ``"pending"`` — mutate a row's fields, then call ``live.update(render())``
    to repaint. *live* has already been ``.start()``-ed.
    """
    from rich.live import Live
    from rich.table import Table

    rows: list[dict[str, str]] = [
        {
            "idx": str(i),
            "label": label,
            "code": step.spec.code,
            "status": "pending",
            "elapsed": "",
            "sites": "",
            "cached": "",
        }
        for i, (label, step) in enumerate(steps, start=1)
    ]

    def render() -> Table:
        table = Table(title=f"{pipeline_name}  ({len(rows)} step(s))")
        table.add_column("#", justify="right")
        table.add_column("Step")
        table.add_column("Code")
        table.add_column("Status")
        table.add_column("Time", justify="right")
        table.add_column("Sites")
        table.add_column("Cached")
        for r in rows:
            style = _STATUS_STYLE.get(r["status"])
            status_text = f"[{style}]{r['status']}[/{style}]" if style else r["status"]
            table.add_row(
                r["idx"], r["label"], r["code"], status_text,
                r["elapsed"], r["sites"], r["cached"],
            )
        return table

    live = Live(render(), refresh_per_second=4, transient=False)
    live.start()
    return live, rows, render
