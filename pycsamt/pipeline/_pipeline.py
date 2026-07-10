"""Core Pipeline class for the pyCSAMT processing engine.

A :class:`Pipeline` is an ordered sequence of :class:`~._steps.Step` objects
that is applied to a :class:`~pycsamt.site.base.Sites` collection in one
call to :meth:`Pipeline.run`.

Design goals
------------
* **scikit-learn-style API** — steps are stored as ``(label, Step)`` tuples;
  the pipeline prints as a formatted table before running.
* **Immutable-while-running** — once :meth:`run` is called the step list is
  frozen until the run completes.
* **Runs to completion** — errors are collected per step; execution continues
  unless ``on_step_error="raise"`` is set in
  :data:`~pycsamt.api.pipe.PYCSAMT_PIPE`.
* **Reproducible** — the canonical YAML config is saved to the output
  directory alongside the results.

Examples
--------
Build and run from code::

    from pycsamt.emtools.pipe import Pipeline, Step

    pipe = Pipeline([
        ("notch",      Step("NR001", mains_hz=50)),
        ("band",       Step("FREQ001")),
        ("align",      Step("FREQ004")),
        ("correct_ss", Step("SS001")),
        ("rotate",     Step("TZ001")),
    ])
    print(pipe)
    result = pipe.run(sites, outdir="willy_results/")

Build from preset and customise::

    pipe = Pipeline.from_preset("full_processing")
    pipe.remove("mask_skew")
    pipe.append("skew_band", Step("SK002", threshold=0.25))

Build from a YAML config::

    pipe = Pipeline.from_yaml("processing/willy_workflow.yaml")
    result = pipe.run(sites)
"""

from __future__ import annotations

import copy
import time
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ._base import PipelineBase
from ._steps import Step, StepResult

# Sentinel: distinguishes "user did not pass outdir" from "user passed None".
_UNSET = object()


# ---------------------------------------------------------------------------
# PipelineResult
# ---------------------------------------------------------------------------

@dataclass
class PipelineResult:
    """Return value of :meth:`Pipeline.run`.

    Attributes
    ----------
    sites_in:
        Original Sites passed to :meth:`Pipeline.run`.
    sites_out:
        Fully processed Sites after all steps.
    step_results:
        One :class:`~._steps.StepResult` per pipeline step, in order.
    outdir:
        Root output directory (``None`` when no output was requested).
    elapsed_sec:
        Total wall-clock time for the run.
    processed_paths:
        Paths of EDI files written to ``<outdir>/processed/``.
    pipeline_name:
        Name label of the pipeline.
    """

    sites_in: Any
    sites_out: Any
    step_results: list[StepResult]
    outdir: Path | None
    elapsed_sec: float
    processed_paths: list[Path] = field(default_factory=list)
    pipeline_name: str = "pipeline"

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def plots(self) -> list[Path]:
        """All figure paths generated across every step."""
        return [p for sr in self.step_results for p in sr.plots]

    @property
    def n_errors(self) -> int:
        """Number of steps that raised an error."""
        return sum(1 for sr in self.step_results if not sr.ok)

    @property
    def ok(self) -> bool:
        """``True`` when every step completed without error."""
        return self.n_errors == 0

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Return a compact text summary of the run."""
        n_in = _count_sites(self.sites_in)
        n_out = _count_sites(self.sites_out)
        lines = [
            f"PipelineResult  '{self.pipeline_name}'",
            f"  Sites   : {n_in} in → {n_out} out",
            f"  Steps   : {len(self.step_results)} "
            f"({len(self.step_results) - self.n_errors} ok, "
            f"{self.n_errors} err)",
            f"  Time    : {self.elapsed_sec:.2f} s",
            f"  Plots   : {len(self.plots)}",
        ]
        if self.outdir is not None:
            lines.append(f"  Output  : {self.outdir}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"PipelineResult(steps={len(self.step_results)}, "
            f"ok={self.ok}, elapsed={self.elapsed_sec:.2f}s, "
            f"plots={len(self.plots)})"
        )


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

class Pipeline(PipelineBase):
    """An ordered, configurable MT processing pipeline.

    Parameters
    ----------
    steps:
        Either a list of ``(label, Step)`` tuples or a plain list of
        :class:`~._steps.Step` objects (labels are auto-generated from the
        step name).

    Examples
    --------
    >>> from pycsamt.emtools.pipe import Pipeline, Step
    >>> pipe = Pipeline([
    ...     ("notch",    Step("NR001")),
    ...     ("band",     Step("FREQ001", band_hz=(0.001, 10000))),
    ...     ("align",    Step("FREQ004")),
    ... ])
    >>> print(pipe)
    """

    def __init__(
        self,
        steps: list[tuple[str, Step]] | list[Step],
        *,
        name: str = "pipeline",
        _output_dir: str | Path | None = None,
    ) -> None:
        self.name = name
        self._output_dir = _output_dir
        self._steps: list[tuple[str, Step]] = self._normalise(steps)
        self._running = False

    # ------------------------------------------------------------------
    # Building the pipeline
    # ------------------------------------------------------------------

    def append(self, label: str, step: Step) -> Pipeline:
        """Add *step* to the end of the pipeline.

        Returns *self* so calls can be chained.
        """
        self._check_mutable("append")
        self._steps.append((label, step))
        return self

    def insert(self, idx: int, label: str, step: Step) -> Pipeline:
        """Insert *step* at position *idx* (0-based).

        Returns *self*.
        """
        self._check_mutable("insert")
        self._steps.insert(idx, (label, step))
        return self

    def remove(self, label: str) -> Pipeline:
        """Remove the first step whose label matches *label*.

        Returns *self*.

        Raises
        ------
        KeyError
            When no step with that label exists.
        """
        self._check_mutable("remove")
        for i, (lbl, _) in enumerate(self._steps):
            if lbl == label:
                del self._steps[i]
                return self
        raise KeyError(f"No step labelled {label!r} in the pipeline.")

    def replace(self, label: str, step: Step) -> Pipeline:
        """Replace the step labelled *label* with *step*.

        Returns *self*.
        """
        self._check_mutable("replace")
        for i, (lbl, _) in enumerate(self._steps):
            if lbl == label:
                self._steps[i] = (label, step)
                return self
        raise KeyError(f"No step labelled {label!r} in the pipeline.")

    def clone(self) -> Pipeline:
        """Return a deep copy of this pipeline."""
        return copy.deepcopy(self)

    # ------------------------------------------------------------------
    # Running
    # ------------------------------------------------------------------

    def run(
        self,
        sites: Any,
        *,
        outdir: Any = _UNSET,
        save_plots: bool = True,
        save_edis: bool = True,
        save_report: bool = True,
        api: Any = None,
    ) -> PipelineResult:
        """Run all steps in order and return a :class:`PipelineResult`.

        Parameters
        ----------
        sites:
            Input :class:`~pycsamt.site.base.Sites` collection.
        outdir:
            Root output directory.  Falls back first to the pipeline's own
            ``_output_dir`` (set by ``from_yaml``), then to
            :attr:`~pycsamt.api.pipe.PipelineAPIConfig.output_root`.
        save_plots:
            Generate and save QC figures for each step.
        save_edis:
            Write processed EDI files after the final step.
        save_report:
            Write HTML and/or text reports to *outdir*.
        api:
            Optional :class:`~pycsamt.api.pipe.PipelineAPIConfig` override.
            When ``None`` the global singleton is used.
        """
        cfg = api or _get_cfg()

        # Resolve the output root:
        #   outdir=_UNSET → use pipeline default → fall back to cfg.output_root
        #   outdir=None   → explicit opt-out, no files written
        #   outdir="path" → use that path
        if outdir is _UNSET:
            root: Any = self._output_dir or cfg.output_root
        else:
            root = outdir  # may be None (opt-out) or an explicit path

        # Resolve output directory
        out = None
        if root is not None:
            from ._output import OutputDir
            out = OutputDir(root, api=cfg)
            out.setup()

        self._running = True
        t_start = time.perf_counter()
        current_sites = sites
        step_results: list[StepResult] = []

        # Save initial pipeline config for reproducibility
        yaml_str = self.to_yaml_string()
        if out is not None:
            out.save_pipeline_config(yaml_str)

        # ── Main loop ─────────────────────────────────────────────────
        for step_idx, (label, step) in enumerate(self._steps, start=1):
            n_in = _count_sites(current_sites)
            t_step = time.perf_counter()
            error: Exception | None = None
            plot_paths: list[Path] = []
            sites_after = current_sites

            # Progress output
            if cfg.show_progress and cfg.progress_style != "silent":
                _print_step_start(step_idx, len(self._steps), label,
                                   step.spec.code)

            # --- Transform -------------------------------------------
            try:
                sites_after = step.transform(current_sites)
            except Exception as exc:
                error = exc
                if cfg.on_step_error == "raise":
                    self._running = False
                    raise
                elif cfg.on_step_error == "warn":
                    warnings.warn(
                        f"Pipeline step {label!r} [{step.spec.code}] "
                        f"raised {type(exc).__name__}: {exc}.  "
                        "Continuing with previous sites.",
                        stacklevel=2,
                    )
                # "skip" — silent, continue with current_sites
                sites_after = current_sites

            # --- QC plots -------------------------------------------
            if save_plots and out is not None and error is None:
                for fn_name, fig in step.generate_qc_plots(sites_after):
                    p = out.save_figure(
                        fig, fn_name, step_idx, label, api=cfg
                    )
                    if p is not None:
                        plot_paths.append(p)
                    try:
                        import matplotlib.pyplot as plt
                        plt.close(fig)
                    except Exception:
                        pass

            # --- Intermediate EDI snapshot --------------------------
            if cfg.save_intermediate and out is not None and error is None:
                snap_dir = out.step_plot_dir(step_idx, label) / "edi_snapshot"
                snap_dir.mkdir(exist_ok=True)
                try:
                    from pycsamt.site.export import (
                        write_sites,
                    )
                    write_sites(sites_after, snap_dir, exist_ok=True)
                except Exception:
                    pass

            elapsed = time.perf_counter() - t_step
            n_out = _count_sites(sites_after)

            sr = StepResult(
                step_idx=step_idx,
                step_name=label,
                step_code=step.spec.code,
                step_label=step.spec.label,
                params=dict(step.params),
                elapsed_sec=elapsed,
                plots=plot_paths,
                n_sites_in=n_in,
                n_sites_out=n_out,
                error=error,
            )
            step_results.append(sr)
            current_sites = sites_after

            if cfg.show_progress and cfg.progress_style != "silent":
                _print_step_done(sr)

        # ── Post-run ──────────────────────────────────────────────────
        total_elapsed = time.perf_counter() - t_start

        # Write final processed EDIs
        processed_paths: list[Path] = []
        if save_edis and out is not None:
            processed_paths = out.write_edis(current_sites)

        # Write reports
        if save_report and out is not None:
            from ._report import (
                make_html_report,
                make_text_report,
            )
            n_in_total = _count_sites(sites)
            n_out_total = _count_sites(current_sites)
            txt = make_text_report(
                self.name, step_results, total_elapsed, out.root,
                n_in_total, n_out_total,
            )
            html = make_html_report(
                self.name, step_results, total_elapsed, out.root,
                n_in_total, n_out_total, pipeline_yaml=yaml_str,
            )
            fmts = cfg.report_formats or ("html", "txt")
            if "txt" in fmts:
                out.save_text(txt, "summary.txt")
            if "html" in fmts:
                out.save_text(html, "report.html")

        self._running = False

        result = PipelineResult(
            sites_in=sites,
            sites_out=current_sites,
            step_results=step_results,
            outdir=out.root if out is not None else None,
            elapsed_sec=total_elapsed,
            processed_paths=processed_paths,
            pipeline_name=self.name,
        )

        if cfg.show_progress and cfg.progress_style != "silent":
            _print_run_done(result)

        return result

    # ------------------------------------------------------------------
    # Config I/O
    # ------------------------------------------------------------------

    @classmethod
    def from_yaml(cls, path: str | Path) -> Pipeline:
        """Load a pipeline from a YAML config file."""
        from ._config import _pipeline_from_dict, load_yaml
        raw = load_yaml(path)
        steps, name, output_dir = _pipeline_from_dict(raw)
        return cls(steps, name=name, _output_dir=output_dir)

    @classmethod
    def from_json(cls, path: str | Path) -> Pipeline:
        """Load a pipeline from a JSON config file."""
        from ._config import _pipeline_from_dict, load_json
        raw = load_json(path)
        steps, name, output_dir = _pipeline_from_dict(raw)
        return cls(steps, name=name, _output_dir=output_dir)

    @classmethod
    def from_py(cls, path: str | Path) -> Pipeline:
        """Load a pipeline from a Python config file.

        The file must expose a ``pipeline_config`` dict.
        """
        from ._config import _pipeline_from_dict, load_py
        raw = load_py(path)
        steps, name, output_dir = _pipeline_from_dict(raw)
        return cls(steps, name=name, _output_dir=output_dir)

    @classmethod
    def from_preset(
        cls,
        name: str,
        pipeline_name: str | None = None,
    ) -> Pipeline:
        """Build a pipeline from a named preset.

        Parameters
        ----------
        name:
            Preset identifier, e.g. ``"full_processing"``.
        pipeline_name:
            Optional override for the pipeline label.

        See also
        --------
        pycsamt.emtools.pipe.preset_catalogue
        """
        from ._presets import get_preset
        preset = get_preset(name)
        return cls(
            list(preset.steps),
            name=pipeline_name or preset.name,
        )

    def to_yaml_string(self) -> str:
        """Serialise this pipeline to a YAML string."""
        from ._config import pipeline_to_yaml
        return pipeline_to_yaml(
            self._steps,
            name=self.name,
            output_dir=str(self._output_dir) if self._output_dir else None,
        )

    def to_yaml(self, path: str | Path) -> None:
        """Write this pipeline config to *path* as YAML."""
        yaml_str = self.to_yaml_string()
        Path(path).write_text(yaml_str, encoding="utf-8")

    def to_json(self, path: str | Path) -> None:
        """Write this pipeline config to *path* as JSON."""
        import json as _json

        import yaml

        from ._config import pipeline_to_yaml

        data_str = pipeline_to_yaml(
            self._steps,
            name=self.name,
            output_dir=str(self._output_dir) if self._output_dir else None,
        )
        try:
            data = yaml.safe_load(data_str)
        except Exception:
            data = {"name": self.name, "steps": [s.to_dict() for _, s in self._steps]}
        Path(path).write_text(
            _json.dumps(data, indent=2, default=str),
            encoding="utf-8",
        )

    # ------------------------------------------------------------------
    # Inspection
    # ------------------------------------------------------------------

    def describe(self) -> Any:
        """Return a :class:`pandas.DataFrame` describing the pipeline steps."""
        try:
            import pandas as pd
        except ImportError:
            return self.__repr__()
        rows = []
        for idx, (label, step) in enumerate(self._steps, start=1):
            rows.append({
                "#": idx,
                "label": label,
                "code": step.spec.code,
                "name": step.spec.name,
                "category": step.spec.category,
                "label_long": step.spec.label,
                "params": step.params,
                "returns_sites": step.spec.returns_sites,
            })
        return pd.DataFrame(rows).set_index("#")

    def steps_in_category(self, category: str) -> list[tuple[str, Step]]:
        """Return steps belonging to *category*."""
        return [(lbl, s) for lbl, s in self._steps
                if s.spec.category == category]

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return _format_repr(self._steps, self.name)

    def _repr_html_(self) -> str:
        return _format_html(self._steps, self.name)

    def __len__(self) -> int:
        return len(self._steps)

    def __iter__(self):
        return iter(self._steps)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._steps[key]
        for lbl, step in self._steps:
            if lbl == key:
                return step
        raise KeyError(key)

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    @staticmethod
    def _normalise(steps: Any) -> list[tuple[str, Step]]:
        """Accept ``[Step, ...]`` or ``[(label, Step), ...]``."""
        normalised = []
        for item in steps:
            if isinstance(item, tuple) and len(item) == 2:
                label, step = item
                if not isinstance(step, Step):
                    raise TypeError(
                        f"Expected (str, Step) tuple, got {type(step)}"
                    )
                normalised.append((str(label), step))
            elif isinstance(item, Step):
                normalised.append((item.spec.name, item))
            else:
                raise TypeError(
                    f"Pipeline steps must be Step instances or (str, Step) "
                    f"tuples.  Got: {type(item)}"
                )
        return normalised

    def _check_mutable(self, op: str) -> None:
        if self._running:
            raise RuntimeError(
                f"Cannot {op!r} a pipeline that is currently running."
            )


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _count_sites(sites: Any) -> int:
    try:
        return len(sites)
    except (TypeError, AttributeError):
        return 0


def _get_cfg() -> Any:
    from pycsamt.api.pipe import PYCSAMT_PIPE
    return PYCSAMT_PIPE


def _print_step_start(idx: int, total: int, label: str, code: str) -> None:
    print(f"  [{idx}/{total}] {label} [{code}] ...", end="", flush=True)


def _print_step_done(sr: StepResult) -> None:
    status = "done" if sr.ok else "ERROR"
    print(f" {status} ({sr.elapsed_sec:.2f}s)", flush=True)


def _print_run_done(result: PipelineResult) -> None:
    status = "OK" if result.ok else f"{result.n_errors} error(s)"
    print(
        f"\nPipeline '{result.pipeline_name}' finished  "
        f"[{status}]  {result.elapsed_sec:.2f}s  "
        f"plots={len(result.plots)}",
        flush=True,
    )


def _format_repr(
    steps: list[tuple[str, Step]],
    name: str,
) -> str:
    """scikit-learn-style multi-line repr."""
    try:
        from pycsamt.api.pipe import PYCSAMT_PIPE
        width = PYCSAMT_PIPE.repr_width
    except Exception:
        width = 80

    n = len(steps)
    title = f"Pipeline  '{name}'"
    header_pad = max(0, width - len(title) - len(f"  {n} steps"))
    header = f"{title}  {'─' * header_pad}  {n} step{'s' if n != 1 else ''}"

    if not steps:
        return header + "\n  (empty)\n" + "─" * width

    # Compute column widths
    max_lbl = max(len(lbl) for lbl, _ in steps)
    max_code = max(len(s.spec.code) for _, s in steps)
    max_label = max(len(s.spec.label) for _, s in steps)

    lbl_w = max(max_lbl, 8)
    code_w = max(max_code + 2, 8)   # +2 for brackets
    label_w = max(max_label, 20)

    lines = [header]
    for idx, (lbl, step) in enumerate(steps, start=1):
        code_col = f"[{step.spec.code}]"
        params_kv = [(k, v) for k, v in step.params.items()]
        params_str = ""
        if params_kv:
            params_str = "  " + "  ".join(f"{k}={v!r}" for k, v in params_kv)
        line = (
            f"  ({idx:2d}) {lbl:<{lbl_w}}  "
            f"{code_col:<{code_w}}  "
            f"{step.spec.label:<{label_w}}"
            f"{params_str}"
        )
        lines.append(line[:width + 20])  # allow slight overflow for params

    lines.append("─" * width)
    return "\n".join(lines)


def _format_html(
    steps: list[tuple[str, Step]],
    name: str,
) -> str:
    """Jupyter-friendly HTML repr."""
    category_colors = {
        "frequency":      "#d4e6f1",
        "noise_removal":  "#d5f5e3",
        "static_shift":   "#fdebd0",
        "tensor":         "#e8daef",
        "dimensionality": "#fdfefe",
        "skew":           "#f9ebea",
        "source_effects": "#eafaf1",
        "qc":             "#f2f3f4",
    }
    rows = ""
    for idx, (lbl, step) in enumerate(steps, start=1):
        bg = category_colors.get(step.spec.category, "#ffffff")
        params_str = (
            ", ".join(f"{k}={v!r}" for k, v in step.params.items())
            if step.params else "—"
        )
        rows += (
            f"<tr style='background:{bg}'>"
            f"<td>{idx}</td>"
            f"<td><b>{lbl}</b></td>"
            f"<td><code>{step.spec.code}</code></td>"
            f"<td>{step.spec.label}</td>"
            f"<td>{step.spec.category}</td>"
            f"<td><code>{params_str}</code></td>"
            "</tr>"
        )
    return (
        f"<b>Pipeline</b> <i>'{name}'</i>  — {len(steps)} step(s)<br/>"
        "<table style='border-collapse:collapse;font-size:0.9em'>"
        "<tr style='background:#4a6fa5;color:white'>"
        "<th>#</th><th>Label</th><th>Code</th>"
        "<th>Description</th><th>Category</th><th>Params</th>"
        "</tr>"
        + rows
        + "</table>"
    )


__all__ = ["Pipeline", "PipelineResult"]
