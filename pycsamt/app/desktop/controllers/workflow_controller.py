# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
WorkflowController — Qt-free model behind the desktop Pipeline window.

It is a thin layer over the library pipeline engine (:mod:`pycsamt.pipeline`),
so the desktop builds and runs exactly what the library, the CLI and saved
YAML/JSON workflows describe:

* **catalogue** — the step registry (55+ steps in categories) and presets;
* **workflow** — an ordered list of :class:`WorkflowStep` (a library
  :class:`~pycsamt.pipeline.Step` + an enabled flag + run state);
* **execution** — :meth:`execute_step` runs one step with the library's own
  ``Step.transform`` and the ``PYCSAMT_PIPE.on_step_error`` policy.  The
  window drives it step by step (live status, Stop between steps, a
  per-step snapshot for previews) instead of calling ``Pipeline.run``,
  whose ``on_step`` hook only fires after a step and cannot stop a run;
* **outputs** — :meth:`build_result` assembles a real
  :class:`~pycsamt.pipeline.PipelineResult`, so the library's reports,
  dashboard plots and run history work unchanged.
"""

from __future__ import annotations

import ast
import inspect
import time
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any


class RunStatus(Enum):
    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    DONE = "done"
    ERROR = "error"
    DISABLED = "disabled"
    OUTDATED = "outdated"


# Filled "pill" colours, all ≥ 4.5 : 1 contrast against white text so they
# read on both the light and dark themes (the old stepper used plain
# green/yellow *text*; yellow was near-invisible on the light theme).
STATUS_STYLE: dict[RunStatus, tuple[str, str]] = {
    RunStatus.PENDING: ("Pending", "#5f6b7a"),
    RunStatus.QUEUED: ("Queued", "#3b5bdb"),
    RunStatus.RUNNING: ("Running", "#1864ab"),
    RunStatus.DONE: ("Done", "#2a7f3f"),  # 5.0:1 (#2b8a3e was 4.37:1)
    RunStatus.ERROR: ("Error", "#c92a2a"),
    RunStatus.DISABLED: ("Off", "#6c757d"),
    RunStatus.OUTDATED: ("Outdated", "#9c4f00"),
}

# Arguments of step functions that are plumbing, not user parameters.
_HIDDEN_ARGS = {
    "sites", "site", "edis", "self", "inplace", "verbose", "recursive",
    "on_dup", "ax", "axes", "fig", "show", "savefig", "figsize", "return_df",
    "return_fig", "copy", "kwargs", "args",
}


@dataclass
class ParamField:
    """One editable parameter of a step."""

    name: str
    value: Any
    default: Any
    kind: str  # "int" | "float" | "bool" | "str" | "literal"
    advanced: bool  # False: registry default; True: from the signature


@dataclass
class WorkflowStep:
    """A library Step plus its desktop run state."""

    step: Any  # pycsamt.pipeline.Step
    label: str
    enabled: bool = True
    status: RunStatus = RunStatus.PENDING
    sites_after: Any = None
    n_in: int = 0
    n_out: int = 0
    elapsed: float = 0.0
    error: str = ""
    result: Any = None  # library StepResult of the last run

    @property
    def code(self) -> str:
        return self.step.spec.code

    @property
    def spec(self):
        return self.step.spec


def _kind_of(value: Any) -> str:
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, str):
        return "str"
    return "literal"


def parse_literal(text: str) -> Any:
    """Parse a free-text parameter (``None``, ``(1e-3, 1e2)``, ``[1, 2]``)."""
    text = text.strip()
    if text == "":
        return None
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return text


def count_sites(sites: Any) -> int:
    if sites is None:
        return 0
    try:
        return len(sites)
    except TypeError:
        try:
            return len(list(sites))
        except Exception:
            return 0


class WorkflowController:
    """Model of the editable, runnable workflow shown by PipelineWindow."""

    def __init__(self) -> None:
        self.steps: list[WorkflowStep] = []
        self.name = "workflow"
        self.input_sites: Any = None
        self.input_label = ""
        self.run_elapsed = 0.0

    # ── Catalogue ─────────────────────────────────────────────────────────────

    @staticmethod
    def categories() -> list[str]:
        from pycsamt.pipeline import categories

        return list(categories())

    @staticmethod
    def catalogue(category: str | None = None) -> list:
        from pycsamt.pipeline import list_steps

        return list(list_steps(category)) if category else list(list_steps())

    @staticmethod
    def presets() -> list:
        from pycsamt.pipeline import list_presets

        return list(list_presets())

    @staticmethod
    def enable_ai_steps() -> list[str]:
        """Register the opt-in AI steps (lazy: no torch import until run)."""
        from pycsamt.pipeline import register_ai_steps

        try:
            return [s.code for s in register_ai_steps()]
        except ValueError:  # already registered
            return []

    @staticmethod
    def describe(spec) -> str:
        """First paragraph of the step function's docstring (or '')."""
        try:
            fn = spec.get_fn()
        except Exception:
            return ""
        doc = inspect.getdoc(fn) or ""
        return doc.split("\n\n", 1)[0].replace("\n", " ").strip()

    # ── Editing ───────────────────────────────────────────────────────────────

    def _unique_label(self, base: str) -> str:
        labels = {s.label for s in self.steps}
        if base not in labels:
            return base
        i = 2
        while f"{base}_{i}" in labels:
            i += 1
        return f"{base}_{i}"

    def add_step(self, code: str, index: int | None = None, **params) -> int:
        from pycsamt.pipeline import Step

        step = Step(code, **params)
        ws = WorkflowStep(step=step, label=self._unique_label(step.spec.name))
        if index is None or index >= len(self.steps):
            self.steps.append(ws)
            index = len(self.steps) - 1
        else:
            self.steps.insert(max(index, 0), ws)
        self._invalidate_from(index)
        return index

    def remove_step(self, index: int) -> None:
        if 0 <= index < len(self.steps):
            del self.steps[index]
            self._invalidate_from(index)

    def move_step(self, index: int, delta: int) -> int:
        new = index + delta
        if not (0 <= index < len(self.steps) and 0 <= new < len(self.steps)):
            return index
        self.steps[index], self.steps[new] = self.steps[new], self.steps[index]
        self._invalidate_from(min(index, new))
        return new

    def set_enabled(self, index: int, enabled: bool) -> None:
        ws = self.steps[index]
        ws.enabled = enabled
        self._invalidate_from(index)

    def clear(self) -> None:
        self.steps = []
        self.run_elapsed = 0.0

    def load_preset(self, name: str) -> None:
        from pycsamt.pipeline import Pipeline

        self.load_pipeline(Pipeline.from_preset(name))

    def load_pipeline(self, pipeline) -> None:
        self.steps = [WorkflowStep(step=st, label=lbl) for lbl, st in pipeline]
        self.name = getattr(pipeline, "name", "workflow") or "workflow"
        self.run_elapsed = 0.0

    def to_pipeline(self, enabled_only: bool = True):
        from pycsamt.pipeline import Pipeline, Step

        steps = [
            (ws.label, Step(ws.code, **ws.step.params))
            for ws in self.steps
            if ws.enabled or not enabled_only
        ]
        return Pipeline(steps, name=self.name)

    def save(self, path: str | Path) -> None:
        """Save as library YAML/JSON (runnable by the CLI and scripts)."""
        path = Path(path)
        pipe = self.to_pipeline(enabled_only=True)
        if path.suffix.lower() == ".json":
            pipe.to_json(path)
        else:
            pipe.to_yaml(path)

    def load(self, path: str | Path) -> None:
        from pycsamt.pipeline import Pipeline

        path = Path(path)
        loader = Pipeline.from_json if path.suffix.lower() == ".json" else (
            Pipeline.from_yaml
        )
        self.load_pipeline(loader(path))

    # ── Parameters ────────────────────────────────────────────────────────────

    def param_fields(self, index: int) -> list[ParamField]:
        """Registry defaults first, then the function's other keywords."""
        ws = self.steps[index]
        spec, params = ws.spec, ws.step.params
        fields_: list[ParamField] = []
        for name, default in spec.defaults.items():
            fields_.append(ParamField(
                name, params.get(name, default), default, _kind_of(default),
                advanced=False,
            ))
        try:
            sig = inspect.signature(spec.get_fn())
        except Exception:
            sig = None
        if sig is not None:
            for name, p in sig.parameters.items():
                if (
                    name in spec.defaults
                    or name in _HIDDEN_ARGS
                    or p.default is inspect.Parameter.empty
                    or p.kind in (p.VAR_KEYWORD, p.VAR_POSITIONAL)
                ):
                    continue
                if not isinstance(p.default, (int, float, str, bool, tuple,
                                              list, type(None))):
                    continue
                fields_.append(ParamField(
                    name, params.get(name, p.default), p.default,
                    _kind_of(p.default), advanced=True,
                ))
        return fields_

    def set_param(self, index: int, name: str, value: Any) -> None:
        ws = self.steps[index]
        default = next(
            (f.default for f in self.param_fields(index) if f.name == name),
            None,
        )
        if value == default and name not in ws.spec.defaults:
            ws.step.params.pop(name, None)  # keep configs minimal
        else:
            ws.step.params[name] = value
        self._invalidate_from(index)

    def reset_params(self, index: int) -> None:
        ws = self.steps[index]
        ws.step.params = dict(ws.spec.defaults)
        self._invalidate_from(index)

    # ── Run state ─────────────────────────────────────────────────────────────

    def set_input(self, sites: Any, label: str = "") -> None:
        self.input_sites = sites
        self.input_label = label
        self._invalidate_from(0)

    def _invalidate_from(self, index: int) -> None:
        """Results at/after *index* no longer match the workflow."""
        for ws in self.steps[max(index, 0):]:
            if not ws.enabled:
                ws.status = RunStatus.DISABLED
            elif ws.status in (RunStatus.DONE, RunStatus.ERROR):
                ws.status = RunStatus.OUTDATED
            elif ws.status is RunStatus.DISABLED:
                ws.status = RunStatus.PENDING

    def reset_run(self) -> None:
        for ws in self.steps:
            ws.status = RunStatus.PENDING if ws.enabled else RunStatus.DISABLED
            ws.sites_after = None
            ws.result = None
            ws.error = ""
            ws.elapsed = 0.0
        self.run_elapsed = 0.0

    def input_for(self, index: int) -> Any:
        """Sites feeding step *index*: the last finished enabled step before
        it, else the workflow input. ``None`` if an earlier step is not done.
        """
        for ws in reversed(self.steps[:index]):
            if not ws.enabled:
                continue
            if ws.status is RunStatus.DONE or (
                ws.status is RunStatus.ERROR and ws.sites_after is not None
            ):
                return ws.sites_after
            return None
        return self.input_sites

    def plan(self, mode: str, index: int = 0) -> list[int]:
        """Step indices to run: ``"all"``, ``"from"`` (index..end), ``"one"``."""
        if mode == "all":
            rng = range(len(self.steps))
        elif mode == "from":
            rng = range(index, len(self.steps))
        else:
            rng = [index]
        return [i for i in rng if 0 <= i < len(self.steps)
                and self.steps[i].enabled]

    def execute_step(self, index: int, sites_in: Any, step_idx: int = 1):
        """Run one step; return ``(sites_out, ok, stop_run)``.

        Honors ``PYCSAMT_PIPE.on_step_error`` (Settings ▸ Pipeline):
        ``raise`` stops the run, ``warn``/``skip`` continue with the input.
        """
        from pycsamt.api.pipe import PYCSAMT_PIPE
        from pycsamt.pipeline import StepResult

        ws = self.steps[index]
        n_in = count_sites(sites_in)
        t0 = time.perf_counter()
        error = None
        try:
            sites_out = ws.step.transform(sites_in)
        except Exception as exc:
            error = exc
            sites_out = sites_in
        ws.elapsed = time.perf_counter() - t0
        ws.n_in, ws.n_out = n_in, count_sites(sites_out)
        ws.sites_after = sites_out
        ws.error = "" if error is None else f"{type(error).__name__}: {error}"
        ws.status = RunStatus.DONE if error is None else RunStatus.ERROR
        ws.result = StepResult(
            step_idx=step_idx, step_name=ws.label, step_code=ws.code,
            step_label=ws.spec.label, params=dict(ws.step.params),
            elapsed_sec=ws.elapsed, n_sites_in=n_in, n_sites_out=ws.n_out,
            error=error,
        )
        stop = error is not None and PYCSAMT_PIPE.on_step_error == "raise"
        return sites_out, error is None, stop

    @property
    def output_sites(self) -> Any:
        for ws in reversed(self.steps):
            if ws.enabled and ws.status in (RunStatus.DONE, RunStatus.ERROR):
                return ws.sites_after
            if ws.enabled:
                return None
        return None

    @property
    def is_complete(self) -> bool:
        enabled = [ws for ws in self.steps if ws.enabled]
        return bool(enabled) and all(
            ws.status in (RunStatus.DONE, RunStatus.ERROR) for ws in enabled
        )

    def build_result(self):
        """A library PipelineResult for reports, dashboards and history."""
        from pycsamt.pipeline import PipelineResult

        results = [ws.result for ws in self.steps
                   if ws.enabled and ws.result is not None]
        return PipelineResult(
            sites_in=self.input_sites,
            sites_out=self.output_sites,
            step_results=results,
            outdir=None,
            elapsed_sec=self.run_elapsed or sum(r.elapsed_sec for r in results),
            pipeline_name=self.name,
        )

    # ── Outputs ───────────────────────────────────────────────────────────────

    def export(self, outdir: str | Path, *, edis: bool = True,
               reports: bool = True, figures: bool = True) -> dict:
        """Write results with the library's own writers (``PYCSAMT_PIPE``
        folders, figure format/DPI and report formats)."""
        from pycsamt.api.pipe import PYCSAMT_PIPE as cfg
        from pycsamt.pipeline._output import OutputDir

        out = OutputDir(outdir, api=cfg)
        out.setup()
        written: dict[str, list] = {"edis": [], "reports": [], "figures": []}
        yaml_str = self.to_pipeline().to_yaml_string()
        out.save_pipeline_config(yaml_str)
        result = self.build_result()

        if figures:
            import matplotlib.pyplot as plt

            enabled = [(k, w) for k, w in enumerate(self.steps) if w.enabled]
            for i, (k, ws) in enumerate(enabled, 1):
                if ws.status is not RunStatus.DONE:
                    continue
                figs = ws.step.generate_qc_plots(
                    ws.sites_after, before=self.input_for(k)
                )
                for fn_name, fig in figs:
                    p = out.save_figure(fig, fn_name, i, ws.label, api=cfg)
                    plt.close(fig)
                    if p is not None:
                        written["figures"].append(p)
                        if ws.result is not None:
                            ws.result.plots.append(p)

        if edis and result.sites_out is not None:
            written["edis"] = list(out.write_edis(result.sites_out))

        if reports:
            from pycsamt.pipeline._report import (
                make_html_report,
                make_text_report,
            )

            args = (self.name, result.step_results, result.elapsed_sec,
                    out.root, count_sites(result.sites_in),
                    count_sites(result.sites_out))
            fmts = cfg.report_formats or ("html", "txt")
            if "txt" in fmts:
                written["reports"].append(
                    out.save_text(make_text_report(*args), "summary.txt"))
            if "html" in fmts:
                written["reports"].append(out.save_text(
                    make_html_report(*args, pipeline_yaml=yaml_str),
                    "report.html"))
            if "dashboard" in fmts:
                from pycsamt.pipeline._dashboard import make_dashboard_html

                written["reports"].append(out.save_text(
                    make_dashboard_html(*args, pipeline_yaml=yaml_str),
                    "dashboard.html"))
        written["reports"] = [p for p in written["reports"] if p]
        written["root"] = out.root
        return written

    def record_history(self, path: str | Path | None = None) -> Path:
        """Append this run to the run-history log (Settings ▸ Pipeline)."""
        from pycsamt.api.pipe import PYCSAMT_PIPE
        from pycsamt.pipeline._history import append_run, default_history_path

        target = Path(path or PYCSAMT_PIPE.history_path or default_history_path())
        append_run(target, self.build_result())
        return target

    @staticmethod
    def load_history(path: str | Path | None = None, last: int = 100) -> list:
        from pycsamt.api.pipe import PYCSAMT_PIPE
        from pycsamt.pipeline import load_history

        target = path or PYCSAMT_PIPE.history_path
        try:
            return list(load_history(target, last=last))
        except Exception:
            return []


__all__ = [
    "ParamField",
    "RunStatus",
    "STATUS_STYLE",
    "WorkflowController",
    "WorkflowStep",
    "count_sites",
    "parse_literal",
]
