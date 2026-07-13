"""Plot helpers for :mod:`pycsamt.pipeline` results.

The pipeline engine stores step timing, status, station counts, and generated
plot paths in :class:`pycsamt.pipeline.PipelineResult`.  This module provides
small Matplotlib helpers for turning that metadata into review figures that
work in notebooks, reports, and the documentation gallery.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ._pipeline import PipelineResult

__all__ = [
    "plot_pipeline_dashboard",
    "plot_pipeline_status",
    "plot_pipeline_timing",
    "plot_site_count_flow",
]

_OK = "#2f9e44"
_ERR = "#c92a2a"
_WARN = "#f08c00"
_BLUE = "#2f6f8f"
_INK = "#243746"
_GRID = "#d7dee5"
_MUTED = "#74808a"


def _require_matplotlib():
    import matplotlib.pyplot as plt

    return plt


def _step_results(result: PipelineResult) -> list[Any]:
    return list(getattr(result, "step_results", []) or [])


def _step_labels(result: PipelineResult) -> list[str]:
    labels = []
    for sr in _step_results(result):
        idx = getattr(sr, "step_idx", len(labels) + 1)
        name = str(getattr(sr, "step_name", f"step_{idx}"))
        labels.append(f"{idx}. {name}")
    return labels


def _short_labels(result: PipelineResult) -> list[str]:
    labels = []
    for sr in _step_results(result):
        idx = getattr(sr, "step_idx", len(labels) + 1)
        name = str(getattr(sr, "step_name", f"step_{idx}"))
        labels.append(f"{idx}\n{name}")
    return labels


def _elapsed(result: PipelineResult) -> np.ndarray:
    vals = []
    for sr in _step_results(result):
        try:
            vals.append(max(float(getattr(sr, "elapsed_sec", 0.0)), 0.0))
        except Exception:
            vals.append(0.0)
    return np.asarray(vals, dtype=float)


def _site_counts(result: PipelineResult) -> tuple[np.ndarray, np.ndarray]:
    n_in, n_out = [], []
    for sr in _step_results(result):
        n_in.append(int(getattr(sr, "n_sites_in", 0) or 0))
        n_out.append(int(getattr(sr, "n_sites_out", 0) or 0))
    return np.asarray(n_in, dtype=float), np.asarray(n_out, dtype=float)


def _plot_counts(result: PipelineResult) -> np.ndarray:
    vals = []
    for sr in _step_results(result):
        plots = getattr(sr, "plots", []) or []
        vals.append(len(plots))
    return np.asarray(vals, dtype=float)


def _status_colors(result: PipelineResult) -> list[str]:
    return [_OK if getattr(sr, "ok", False) else _ERR for sr in _step_results(result)]


def _style_axis(ax: Any, *, xgrid: bool = False, ygrid: bool = True) -> None:
    ax.set_axisbelow(True)
    if ygrid:
        ax.grid(axis="y", color=_GRID, lw=0.7, alpha=0.75)
    if xgrid:
        ax.grid(axis="x", color=_GRID, lw=0.7, alpha=0.55)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#b8c1ca")
    ax.tick_params(colors=_INK, labelsize=8)
    ax.xaxis.label.set_color(_INK)
    ax.yaxis.label.set_color(_INK)
    ax.title.set_color(_INK)


def _empty_axis(ax: Any, title: str) -> Any:
    ax.text(
        0.5,
        0.5,
        "no pipeline steps",
        ha="center",
        va="center",
        transform=ax.transAxes,
        color=_MUTED,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    for spine in ax.spines.values():
        spine.set_color("#d0d7de")
    return ax


def _annotate_bars(ax: Any, bars: Any, values: np.ndarray, fmt: str) -> None:
    finite = values[np.isfinite(values)]
    ymax = float(np.nanmax(finite)) if finite.size else 1.0
    pad = max(ymax * 0.025, 0.02)
    for bar, val in zip(bars, values):
        if not np.isfinite(val):
            continue
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + pad,
            fmt.format(val),
            ha="center",
            va="bottom",
            fontsize=8,
            color=_INK,
        )


def plot_pipeline_status(
    result: PipelineResult,
    *,
    ax: Any | None = None,
    ok_color: str = _OK,
    error_color: str = _ERR,
    annotate_errors: bool = True,
) -> Any:
    """Plot per-step success status for a pipeline run.

    Failed steps are coloured red and, when available, the exception message
    is written below the corresponding bar.  Empty results draw a neutral
    placeholder instead of failing.
    """
    plt = _require_matplotlib()
    if ax is None:
        _, ax = plt.subplots(figsize=(9.5, 4.0))

    labels = _step_labels(result)
    if not labels:
        return _empty_axis(ax, f"Pipeline status: {result.pipeline_name}")

    x = np.arange(len(labels))
    colors = [
        ok_color if getattr(sr, "ok", False) else error_color
        for sr in _step_results(result)
    ]
    bars = ax.bar(x, np.ones(len(labels)), color=colors, edgecolor="#1f2933")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    ax.set_yticks([])
    ax.set_ylim(-0.32 if annotate_errors else 0.0, 1.25)
    ax.set_ylabel("status")
    ax.set_title(f"Pipeline status: {result.pipeline_name}")

    for i, (bar, sr) in enumerate(zip(bars, _step_results(result))):
        ok = bool(getattr(sr, "ok", False))
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            0.52,
            "OK" if ok else "ERR",
            ha="center",
            va="center",
            color="white",
            fontsize=8,
            fontweight="bold",
        )
        err = getattr(sr, "error", None)
        if annotate_errors and err is not None:
            msg = str(err).strip().replace("\n", " ")
            if len(msg) > 34:
                msg = msg[:31] + "..."
            ax.text(
                i,
                -0.08,
                msg or type(err).__name__,
                ha="center",
                va="top",
                rotation=35,
                fontsize=7,
                color=error_color,
            )
    _style_axis(ax, ygrid=False)
    return ax


def plot_pipeline_timing(
    result: PipelineResult,
    *,
    ax: Any | None = None,
    color: str = _BLUE,
    slow_color: str = _WARN,
    slow_quantile: float = 0.80,
    annotate: bool = True,
) -> Any:
    """Plot elapsed time per pipeline step.

    Steps above *slow_quantile* are highlighted, helping users quickly spot
    expensive processing or plotting stages.
    """
    plt = _require_matplotlib()
    if ax is None:
        _, ax = plt.subplots(figsize=(9.5, 4.0))

    labels = _step_labels(result)
    if not labels:
        return _empty_axis(ax, f"Pipeline timing: {result.pipeline_name}")

    elapsed = _elapsed(result)
    finite = elapsed[np.isfinite(elapsed)]
    q = float(np.nanquantile(finite, slow_quantile)) if finite.size else np.inf
    colors = [slow_color if v >= q and v > 0 else color for v in elapsed]
    x = np.arange(len(labels))
    bars = ax.bar(x, elapsed, color=colors, edgecolor="#1f2933")
    if annotate:
        _annotate_bars(ax, bars, elapsed, "{:.2f}s")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Elapsed time (s)")
    ax.set_title(f"Pipeline timing: {result.pipeline_name}")
    ax.set_ylim(0.0, max(float(np.nanmax(elapsed)) * 1.18, 1.0))
    _style_axis(ax)
    return ax


def plot_site_count_flow(
    result: PipelineResult,
    *,
    ax: Any | None = None,
    color_in: str = "#687582",
    color_out: str = "#7c4d79",
    drop_color: str = _ERR,
    annotate: bool = True,
) -> Any:
    """Plot station counts entering and leaving each pipeline step.

    A faint red band marks steps where the output site count is lower than
    the input site count, which is useful for QC workflows that reject stations.
    """
    plt = _require_matplotlib()
    if ax is None:
        _, ax = plt.subplots(figsize=(9.5, 4.0))

    labels = _step_labels(result)
    if not labels:
        return _empty_axis(ax, f"Pipeline station flow: {result.pipeline_name}")

    x = np.arange(len(labels))
    n_in, n_out = _site_counts(result)
    for xi, a, b in zip(x, n_in, n_out):
        if np.isfinite(a) and np.isfinite(b) and b < a:
            ax.axvspan(xi - 0.45, xi + 0.45, color=drop_color, alpha=0.08)
    ax.plot(x, n_in, marker="o", color=color_in, lw=1.8, label="input")
    ax.plot(x, n_out, marker="s", color=color_out, lw=1.8, label="output")
    if annotate:
        for xi, a, b in zip(x, n_in, n_out):
            ax.text(xi, a, f"{int(a)}", ha="center", va="bottom", fontsize=7)
            ax.text(xi, b, f"{int(b)}", ha="center", va="top", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Stations")
    ax.set_title(f"Pipeline station flow: {result.pipeline_name}")
    ax.legend(frameon=False, ncol=2, fontsize=8)
    ymax = max(float(np.nanmax([n_in, n_out])), 1.0)
    ax.set_ylim(-0.05 * ymax, 1.15 * ymax)
    _style_axis(ax)
    return ax


def _plot_output_counts(result: PipelineResult, *, ax: Any) -> Any:
    labels = _step_labels(result)
    if not labels:
        return _empty_axis(ax, "Generated figures")
    plots = _plot_counts(result)
    x = np.arange(len(labels))
    bars = ax.bar(x, plots, color="#4c6ef5", edgecolor="#1f2933")
    _annotate_bars(ax, bars, plots, "{:.0f}")
    ax.set_xticks(x)
    ax.set_xticklabels(_short_labels(result), fontsize=7)
    ax.set_ylabel("Plots")
    ax.set_title("Generated figures by step")
    ax.set_ylim(0.0, max(float(np.nanmax(plots)) * 1.25, 1.0))
    _style_axis(ax)
    return ax


def _plot_summary_cards(result: PipelineResult, *, ax: Any) -> Any:
    ax.axis("off")
    steps = _step_results(result)
    n_steps = len(steps)
    n_ok = sum(1 for sr in steps if getattr(sr, "ok", False))
    n_err = n_steps - n_ok
    elapsed = float(getattr(result, "elapsed_sec", 0.0) or 0.0)
    plots = len(getattr(result, "plots", []) or [])
    paths = len(getattr(result, "processed_paths", []) or [])
    n_in, n_out = _site_counts(result)
    in0 = int(n_in[0]) if n_in.size else 0
    outn = int(n_out[-1]) if n_out.size else 0
    cards = [
        ("Steps", f"{n_ok}/{n_steps} ok", _OK if n_err == 0 else _WARN),
        ("Errors", str(n_err), _OK if n_err == 0 else _ERR),
        ("Sites", f"{in0} -> {outn}", _BLUE),
        ("Time", f"{elapsed:.2f}s", _BLUE),
        ("Figures", str(plots), "#4c6ef5"),
        ("EDI files", str(paths), "#7c4d79"),
    ]
    for i, (title, value, color) in enumerate(cards):
        x0 = (i % 3) / 3.0 + 0.02
        y0 = 0.58 if i < 3 else 0.12
        width = 0.29
        height = 0.30
        rect = plt_rectangle(ax, x0, y0, width, height, color)
        ax.add_patch(rect)
        ax.text(
            x0 + 0.025,
            y0 + height - 0.07,
            title,
            transform=ax.transAxes,
            fontsize=8,
            color=_MUTED,
            va="top",
        )
        ax.text(
            x0 + 0.025,
            y0 + 0.08,
            value,
            transform=ax.transAxes,
            fontsize=14,
            color=_INK,
            fontweight="bold",
            va="bottom",
        )
    ax.set_title(f"Pipeline summary: {result.pipeline_name}", loc="left")
    return ax


def plt_rectangle(ax: Any, x: float, y: float, w: float, h: float, color: str) -> Any:
    from matplotlib.patches import FancyBboxPatch

    return FancyBboxPatch(
        (x, y),
        w,
        h,
        transform=ax.transAxes,
        boxstyle="round,pad=0.008,rounding_size=0.018",
        facecolor="#f8fafc",
        edgecolor=color,
        linewidth=1.4,
    )


def plot_pipeline_dashboard(
    result: PipelineResult,
    *,
    figsize: tuple[float, float] = (12.0, 8.0),
) -> Any:
    """Create a compact dashboard for a completed pipeline run.

    The dashboard combines run-level cards, per-step status, timing,
    station-count flow, and generated-figure counts.  It is the most useful
    single figure to place in a processing report or notebook.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the dashboard.
    """
    plt = _require_matplotlib()
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=(0.85, 1.25, 1.2))
    ax_cards = fig.add_subplot(gs[0, :])
    ax_status = fig.add_subplot(gs[1, 0])
    ax_timing = fig.add_subplot(gs[1, 1])
    ax_sites = fig.add_subplot(gs[2, 0])
    ax_plots = fig.add_subplot(gs[2, 1])

    _plot_summary_cards(result, ax=ax_cards)
    plot_pipeline_status(result, ax=ax_status, annotate_errors=False)
    plot_pipeline_timing(result, ax=ax_timing, annotate=False)
    plot_site_count_flow(result, ax=ax_sites, annotate=False)
    _plot_output_counts(result, ax=ax_plots)
    return fig
