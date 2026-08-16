"""Global pipeline configuration for pyCSAMT.

Centralises all output, display, and execution settings used by
:mod:`pycsamt.pipeline`.  Follows the same singleton + configure/reset
pattern as the rest of :mod:`pycsamt.api`.

Examples
--------
Override output format temporarily::

    from pycsamt.api.pipe import PYCSAMT_PIPE

    with PYCSAMT_PIPE.context(plot_dpi=300, plot_fmt="pdf"):
        result = pipeline.run(sites)

Change global defaults::

    from pycsamt.api.pipe import configure_pipe

    configure_pipe(output_root="my_project/results", show_progress=False)
"""

from __future__ import annotations

import copy
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass, field, fields
from typing import Any

_ON_STEP_ERROR = {"raise", "warn", "skip"}
_PROGRESS_STYLES = {"bar", "log", "silent", "rich"}
_REPORT_FORMATS = {"html", "txt", "dashboard"}


@dataclass
class PipelineAPIConfig:
    """Runtime configuration for the pyCSAMT pipeline engine.

    Attributes
    ----------
    output_root:
        Default root directory when no *outdir* is passed to
        :meth:`~pycsamt.pipeline.Pipeline.run`.
    processed_subdir:
        Sub-directory inside *output_root* that receives the processed EDI
        files.
    plots_subdir:
        Sub-directory inside *output_root* that receives all QC figures,
        organised into one sub-folder per pipeline step.
    on_step_error:
        Behaviour when a step raises an exception.  ``"raise"`` re-raises
        immediately; ``"warn"`` logs a warning and continues with the
        previous sites; ``"skip"`` is silent.
    save_intermediate:
        If ``True``, EDI files are also written after each step (into
        ``plots_subdir/<step>/edi_snapshot/``).
    show_progress:
        Enable console progress output while the pipeline runs.
    progress_style:
        ``"bar"``  – tqdm-style progress bar (requires *tqdm*; falls back to
        ``"log"`` when unavailable).
        ``"log"``  – one line per step.
        ``"silent"`` – no output.
        ``"rich"``  – a live-updating :class:`rich.table.Table` showing every
        step's status (pending/running/OK/ERR/cached) at once, rewritten in
        place as the run progresses.  See :doc:`/user_guide/pipeline/observability`.
    repr_width:
        Character width used when printing the pipeline summary.
    plot_dpi:
        DPI for saved figures.
    plot_fmt:
        File extension / format for saved figures (``"png"``, ``"pdf"``,
        ``"svg"``).
    report_formats:
        Sequence of report types to write.  Each item must be ``"html"``,
        ``"txt"``, or ``"dashboard"``.  ``"dashboard"`` writes a richer,
        branded ``dashboard.html`` with KPI stat tiles and inline-SVG
        charts (step status, per-step duration, site-count flow) — see
        :doc:`/user_guide/pipeline/outputs`.
    cache_root:
        Default root directory for the pipeline step cache
        (:class:`~pycsamt.pipeline.StepCache`) when
        :meth:`~pycsamt.pipeline.Pipeline.run` is called with
        ``cache=True``.  ``None`` resolves lazily to
        ``~/.pycsamt/pipeline_cache``.
    history_path:
        Default file for the run history log
        (:func:`~pycsamt.pipeline.load_history`) when
        :meth:`~pycsamt.pipeline.Pipeline.run` is called with
        ``history=True``.  ``None`` resolves lazily to
        ``~/.pycsamt/pipeline_history.jsonl``.
    """

    output_root: str = "pipe_results"
    processed_subdir: str = "processed"
    plots_subdir: str = "plots"
    on_step_error: str = "warn"
    save_intermediate: bool = False
    show_progress: bool = True
    progress_style: str = "bar"
    repr_width: int = 80
    plot_dpi: int = 150
    plot_fmt: str = "png"
    report_formats: tuple = field(default_factory=lambda: ("html", "txt"))
    cache_root: str | None = None
    history_path: str | None = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def configure(self, **kw: Any) -> None:
        """Set one or more configuration attributes by keyword."""
        for k, v in kw.items():
            if not hasattr(self, k):
                valid = [f.name for f in fields(self)]
                raise AttributeError(
                    f"Unknown pipe config key {k!r}.  Valid keys: {valid}"
                )
            setattr(self, k, v)

    @contextmanager
    def context(self, **kw: Any) -> Generator[PipelineAPIConfig, None, None]:
        """Temporarily override config values, then restore them."""
        snapshot = {k: getattr(self, k) for k in kw}
        try:
            self.configure(**kw)
            yield self
        finally:
            for k, v in snapshot.items():
                setattr(self, k, v)

    def reset(self) -> None:
        """Restore all attributes to package defaults."""
        _defaults = PipelineAPIConfig()
        for f in fields(self):
            setattr(self, f.name, getattr(_defaults, f.name))

    def clone(self) -> PipelineAPIConfig:
        """Return a deep copy of this config."""
        return copy.deepcopy(self)

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        lines = ["PipelineAPIConfig"]
        for f in fields(self):
            lines.append(f"  {f.name}: {getattr(self, f.name)!r}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Singleton + convenience helpers
# ---------------------------------------------------------------------------

PYCSAMT_PIPE: PipelineAPIConfig = PipelineAPIConfig()


def configure_pipe(**kw: Any) -> None:
    """Configure the global :data:`PYCSAMT_PIPE` singleton."""
    PYCSAMT_PIPE.configure(**kw)


def reset_pipe() -> None:
    """Reset :data:`PYCSAMT_PIPE` to package defaults."""
    PYCSAMT_PIPE.reset()


__all__ = [
    "PYCSAMT_PIPE",
    "PipelineAPIConfig",
    "configure_pipe",
    "reset_pipe",
]
