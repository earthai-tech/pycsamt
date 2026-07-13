# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AgentWorker — QThread worker for pycsamt agent / processing execution.

Phase 4: Supports two execution modes:

  * **LLM mode** (type='llm')
      Instantiates the named agent class from ``pycsamt.agents``, attaches
      a logging handler that streams lines back via ``log_line``, and calls
      ``agent.execute({'sites': sites, ...})``.

  * **Processing mode** (type='processing')
      Calls the named emtools function directly (``et.<fn_name>(sites, **params)``).
      Progress is approximate (0 → 50 → 100).

Signals
-------
log_line(str)       Each captured log / print line from the agent
progress(int)       0–100 progress estimate
result(object)      AgentResult (LLM) or raw function return value
error(str)          Human-readable error on failure
"""

from __future__ import annotations

import io
import logging
from contextlib import redirect_stderr, redirect_stdout
from typing import Any

from PySide6.QtCore import QThread, Signal

from pycsamt.app.desktop.agent_registry import get_entry

# ──────────────────────────────────────────────────────────────────────────────
# Logging bridge
# ──────────────────────────────────────────────────────────────────────────────


class _SignalLogHandler(logging.Handler):
    """Forwards log records to a Qt Signal callable (thread-safe)."""

    def __init__(self, emit_fn):
        super().__init__()
        self._emit_fn = emit_fn
        self.setFormatter(logging.Formatter("%(levelname)s  %(message)s"))

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self._emit_fn(self.format(record))
        except Exception:
            pass


class _StreamCapture(io.StringIO):
    """StringIO that also fires a callback on each write."""

    def __init__(self, callback):
        super().__init__()
        self._cb = callback

    def write(self, s: str) -> int:
        for line in s.splitlines():
            line = line.strip()
            if line:
                self._cb(line)
        return len(s)

    def flush(self) -> None:
        pass


# ──────────────────────────────────────────────────────────────────────────────
# Worker thread
# ──────────────────────────────────────────────────────────────────────────────


class AgentWorker(QThread):
    """Background thread for LLM agents and processing operations."""

    log_line = Signal(str)
    progress = Signal(int)
    result = Signal(object)
    error = Signal(str)
    agent_ready = Signal(object)  # emitted with agent instance after creation

    def __init__(
        self,
        agent_name: str,
        sites,
        params: dict[str, Any],
        api_key: str | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._agent_name = agent_name
        self._sites = sites
        self._params = params
        self._api_key = api_key
        self._cancelled = False

    # ── Cancellation ──────────────────────────────────────────────────

    def cancel(self) -> None:
        self._cancelled = True

    # ── Thread entry ──────────────────────────────────────────────────

    def run(self) -> None:
        entry = get_entry(self._agent_name)
        if entry is None:
            self.error.emit(f"Unknown agent: '{self._agent_name}'")
            return

        # Attach logging + stdout capture
        log_handler = _SignalLogHandler(self.log_line.emit)
        root_logger = logging.getLogger()
        root_logger.addHandler(log_handler)
        cap = _StreamCapture(self.log_line.emit)

        try:
            # Suppress plt.show() calls inside emtools/agent functions — we
            # render figures ourselves via MplCanvas.show_figure().
            import unittest.mock as _mock

            import matplotlib.pyplot as _plt

            with (
                redirect_stdout(cap),
                redirect_stderr(cap),
                _mock.patch.object(_plt, "show", lambda *a, **k: None),
            ):
                if entry["type"] == "llm":
                    res = self._run_llm(entry)
                else:
                    res = self._run_processing(entry)
            if not self._cancelled:
                self.result.emit(res)
        except Exception as exc:
            self.error.emit(str(exc))
        finally:
            root_logger.removeHandler(log_handler)

    # ── LLM mode ──────────────────────────────────────────────────────

    def _run_llm(self, entry: dict):
        import pycsamt.agents as ag

        self.log_line.emit(f"Starting {self._agent_name} (LLM mode)…")
        self.progress.emit(5)

        cls = getattr(ag, entry["class_name"])
        init_kwargs = {
            k: v
            for k, v in self._params.items()
            if k in cls.__init__.__code__.co_varnames
        }
        if self._api_key:
            init_kwargs["api_key"] = self._api_key

        agent = cls(**init_kwargs)
        self.agent_ready.emit(agent)

        self.log_line.emit("Running agent.execute()…")
        self.progress.emit(20)

        input_data: dict = {"sites": self._sites}
        # Forward extra params (user_prompt, context, etc.) that are not init
        # kwargs — agents read them from input_data via input_data.get()
        for k, v in self._params.items():
            if k not in init_kwargs and v not in (None, "", False):
                input_data[k] = v
        res = agent.execute(input_data)

        self.progress.emit(100)
        self.log_line.emit(
            f"Done — status: {getattr(res, 'status', '?')}  "
            f"elapsed: {getattr(res, 'elapsed_seconds', 0):.1f}s"
        )
        if getattr(res, "warnings", None):
            for w in res.warnings:
                self.log_line.emit(f"WARNING  {w}")
        return res

    # ── Processing mode ───────────────────────────────────────────────

    def _run_processing(self, entry: dict):
        import inspect

        import pycsamt.emtools as et

        fn_name = entry["fn_name"]
        self.log_line.emit(f"Running et.{fn_name}()…")
        self.progress.emit(10)

        fn = getattr(et, fn_name)

        def _call(func, kwargs_src: dict, data=None):
            """Call func(sites, **filtered_kwargs), adding verbose=0 only if accepted.

            data overrides self._sites as the first positional argument; used when
            a processing step returns corrected Sites that the plot step should use.
            """
            sig_params = set(inspect.signature(func).parameters.keys())
            kw = {k: v for k, v in kwargs_src.items() if k in sig_params}
            if "verbose" in sig_params:
                kw["verbose"] = 0
            sites_arg = data if data is not None else self._sites
            return func(sites_arg, **kw)

        self.progress.emit(30)
        res = _call(fn, self._params)
        self.progress.emit(70)
        self.log_line.emit(f"et.{fn_name}() completed.")

        # If result is not renderable (no figure/axes) AND a separate result_plot
        # function is defined, call it now to produce the figure.
        result_plot_fn = entry.get("result_plot")
        if (
            result_plot_fn
            and result_plot_fn != fn_name
            and not _has_figure(res)
        ):
            plot_fn = getattr(et, result_plot_fn, None)
            if plot_fn is not None:
                try:
                    self.log_line.emit(f"Plotting et.{result_plot_fn}()…")
                    # When the processing step returned corrected Sites, pass it
                    # as the data source for the result_plot (not the original sites).
                    plot_data = res if _is_sites_like(res) else None
                    ax = _call(plot_fn, self._params, data=plot_data)
                    if ax is not None:
                        res = ax  # Axes → panel extracts .get_figure()
                except Exception as exc:
                    self.log_line.emit(f"Plot step failed: {exc}")

        self.progress.emit(100)
        return res


def _has_figure(obj) -> bool:
    """Return True if obj is or contains a matplotlib Figure / Axes."""
    if obj is None:
        return False
    if hasattr(obj, "savefig"):  # Figure
        return True
    if hasattr(obj, "get_figure"):  # Axes
        return True
    if hasattr(obj, "data") and isinstance(obj.data, dict):
        return any(
            hasattr(v, "savefig") or hasattr(v, "get_figure")
            for v in obj.data.values()
        )
    return False


# ──────────────────────────────────────────────────────────────────────────────
# Follow-up chat worker
# ──────────────────────────────────────────────────────────────────────────────


class ChatWorker(QThread):
    """Background thread for post-run follow-up LLM queries.

    Calls ``agent.query_llm()`` with an optional result-context prefix so the
    model can answer questions about the just-completed analysis.
    """

    reply_done = Signal(str)
    error = Signal(str)

    def __init__(self, agent, question: str, context: str = "", parent=None):
        super().__init__(parent)
        self._agent = agent
        self._question = question.strip()
        self._context = context.strip()

    def run(self) -> None:
        try:
            if not getattr(self._agent, "api_key", None):
                self.error.emit("No API key configured for this agent.")
                return
            prompt = self._question
            if self._context:
                prompt = (
                    f"Previous analysis results:\n{self._context}\n\n"
                    f"Follow-up question: {self._question}"
                )
            reply = self._agent.query_llm(prompt) or "(no response from LLM)"
            self.reply_done.emit(reply)
        except Exception as exc:
            self.error.emit(str(exc))


def _is_sites_like(obj) -> bool:
    """Return True if obj looks like a Sites collection."""
    if obj is None:
        return False
    if _has_figure(obj):
        return False
    try:
        from pycsamt.site.base import Sites

        if isinstance(obj, Sites):
            return True
    except Exception:
        pass
    # Structural fallback: Sites exposes edic, by_index, and as_list
    return (
        hasattr(obj, "edic")
        and hasattr(obj, "by_index")
        and hasattr(obj, "as_list")
    )
