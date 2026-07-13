"""Tests for :mod:`pycsamt.pipeline.plot`."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from pycsamt.pipeline import (  # noqa: E402
    PipelineResult,
    StepResult,
    plot_pipeline_dashboard,
    plot_pipeline_status,
    plot_pipeline_timing,
    plot_site_count_flow,
)


def _result(*, failed: bool = False, empty: bool = False) -> PipelineResult:
    if empty:
        steps = []
    else:
        steps = [
            StepResult(
                step_idx=1,
                step_name="notch",
                step_code="NR001",
                step_label="Power-line Harmonic Notch",
                params={},
                elapsed_sec=0.12,
                plots=[Path("notch.png")],
                n_sites_in=5,
                n_sites_out=5,
            ),
            StepResult(
                step_idx=2,
                step_name="mask",
                step_code="QC001",
                step_label="Quality Mask",
                params={},
                elapsed_sec=0.34,
                plots=[],
                n_sites_in=5,
                n_sites_out=4,
                error=RuntimeError("bad station") if failed else None,
            ),
        ]
    return PipelineResult(
        sites_in=None,
        sites_out=None,
        step_results=steps,
        outdir=None,
        elapsed_sec=0.46,
        processed_paths=[Path("a.edi")],
        pipeline_name="demo",
    )


def test_individual_pipeline_plots_return_axes():
    result = _result(failed=True)
    assert plot_pipeline_status(result).get_title()
    assert plot_pipeline_timing(result).get_title()
    assert plot_site_count_flow(result).get_title()


def test_pipeline_dashboard_returns_figure():
    fig = plot_pipeline_dashboard(_result(failed=True))
    assert len(fig.axes) == 5


def test_pipeline_plots_handle_empty_result():
    result = _result(empty=True)
    assert plot_pipeline_status(result).texts
    assert plot_pipeline_timing(result).texts
    assert plot_site_count_flow(result).texts
    assert len(plot_pipeline_dashboard(result).axes) == 5
