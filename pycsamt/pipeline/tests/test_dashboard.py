"""Tests for the "dashboard" report tier (:mod:`pycsamt.pipeline._dashboard`).

Unit tests exercise :func:`make_dashboard_html` directly against small
:class:`~pycsamt.pipeline.StepResult` doubles (built with ``_make_step_result``,
self-contained in this file per this test suite's convention). A separate
integration test at the bottom runs a real pipeline against the WILLY AMT
dataset with ``report_formats=("html", "txt", "dashboard")`` and checks the
actual ``dashboard.html`` written to disk.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from pycsamt.pipeline import StepResult
from pycsamt.pipeline._dashboard import make_dashboard_html

# ─────────────────────────────────────────────────────────────────────────────
# Test doubles
# ─────────────────────────────────────────────────────────────────────────────


def _make_step_result(
    idx: int,
    *,
    ok: bool = True,
    cached: bool = False,
    elapsed_sec: float | None = None,
    n_plots: int = 0,
) -> StepResult:
    return StepResult(
        step_idx=idx,
        step_name=f"step_{idx}",
        step_code=f"NR00{idx}",
        step_label=f"Step {idx} Label",
        params={"k": idx},
        elapsed_sec=elapsed_sec if elapsed_sec is not None else 0.1 * idx,
        plots=[Path(f"plot_{idx}_{j}.png") for j in range(n_plots)],
        n_sites_in=5,
        n_sites_out=5 if ok else 0,
        error=None if ok else RuntimeError("test error"),
        cached=cached,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Structure / branding
# ─────────────────────────────────────────────────────────────────────────────


class TestDashboardStructure:
    def test_returns_valid_html_document(self):
        sr = [_make_step_result(i) for i in range(1, 4)]
        html = make_dashboard_html("my_pipe", sr, 1.23, None, 5, 5)
        assert html.startswith("<!DOCTYPE html>")
        assert "<html" in html
        assert "</html>" in html

    def test_contains_pipeline_name(self):
        sr = [_make_step_result(1)]
        html = make_dashboard_html("willy_pipe", sr, 0.5, None, 5, 5)
        assert "willy_pipe" in html

    def test_contains_step_codes(self):
        sr = [_make_step_result(i) for i in range(1, 3)]
        html = make_dashboard_html("pipe", sr, 0.5, None, 5, 5)
        assert "NR001" in html
        assert "NR002" in html

    def test_success_badge_when_all_ok(self):
        sr = [_make_step_result(1), _make_step_result(2)]
        html = make_dashboard_html("pipe", sr, 0.5, None, 5, 5)
        assert 'status-badge ok">SUCCESS' in html

    def test_error_badge_when_a_step_failed(self):
        sr = [_make_step_result(1, ok=False)]
        html = make_dashboard_html("pipe", sr, 0.5, None, 5, 0)
        assert 'status-badge err">1 ERROR' in html

    def test_dark_mode_media_query_present(self):
        html = make_dashboard_html("pipe", [_make_step_result(1)], 0.1, None, 5, 5)
        assert "@media (prefers-color-scheme: dark)" in html

    def test_scoped_under_dashboard_root_class(self):
        html = make_dashboard_html("pipe", [_make_step_result(1)], 0.1, None, 5, 5)
        assert "pycsamt-dashboard" in html


class TestDashboardBranding:
    def test_real_brand_hex_values_present(self):
        html = make_dashboard_html("pipe", [_make_step_result(1)], 0.1, None, 5, 5)
        # pyCSAMT brand tokens copied from docs/source/_static/css/custom.css
        for hex_value in ("#f15a29", "#3e65b0", "#24324b"):
            assert hex_value in html

    def test_status_palette_hex_values_present(self):
        html = make_dashboard_html("pipe", [_make_step_result(1)], 0.1, None, 5, 5)
        for hex_value in ("#0ca30c", "#fab219", "#d03b3b"):
            assert hex_value in html

    def test_logo_svg_embedded(self):
        html = make_dashboard_html("pipe", [_make_step_result(1)], 0.1, None, 5, 5)
        assert '<svg class="logo"' in html
        assert 'aria-label="pyCSAMT"' in html

    def test_logo_is_the_real_brand_mark_not_a_placeholder(self):
        # The real pycsamt-v2-symbol.svg has this exact viewBox and uses
        # many numbered fill classes (.cls-1..cls-42) for its shading —
        # a stand-in monogram would not reproduce either.
        html = make_dashboard_html("pipe", [_make_step_result(1)], 0.1, None, 5, 5)
        assert 'viewBox="0 0 476.36 454.92"' in html
        assert ".cls-42{fill:none;}" in html

    def test_favicon_link_present(self):
        html = make_dashboard_html("pipe", [_make_step_result(1)], 0.1, None, 5, 5)
        assert '<link rel="icon" type="image/x-icon"' in html
        assert "data:image/x-icon;base64," in html

    def test_icons_use_currentcolor_not_hardcoded_black(self):
        html = make_dashboard_html("pipe", [_make_step_result(1)], 0.1, None, 5, 5)
        assert 'fill="currentColor"' in html
        assert 'fill="#000000"' not in html
        assert 'fill="#1C274C"' not in html


# ─────────────────────────────────────────────────────────────────────────────
# Stat tiles
# ─────────────────────────────────────────────────────────────────────────────


class TestStatTiles:
    def test_shows_ok_over_total_steps(self):
        sr = [_make_step_result(1), _make_step_result(2, ok=False)]
        html = make_dashboard_html("pipe", sr, 1.0, None, 5, 5)
        assert "1/2 ok" in html

    def test_shows_total_time(self):
        html = make_dashboard_html("pipe", [_make_step_result(1)], 4.56, None, 5, 5)
        assert "4.56s" in html

    def test_shows_site_flow_in_tile(self):
        html = make_dashboard_html("pipe", [_make_step_result(1)], 0.1, None, 8, 6)
        assert "8→6" in html

    def test_shows_cache_hit_rate(self):
        sr = [_make_step_result(1, cached=True), _make_step_result(2, cached=False)]
        html = make_dashboard_html("pipe", sr, 1.0, None, 5, 5)
        assert "50%" in html

    def test_figures_tile_counts_all_plots(self):
        sr = [_make_step_result(1, n_plots=2), _make_step_result(2, n_plots=3)]
        html = make_dashboard_html("pipe", sr, 1.0, None, 5, 5)
        assert '<div class="value">5</div>' in html


# ─────────────────────────────────────────────────────────────────────────────
# Charts
# ─────────────────────────────────────────────────────────────────────────────


class TestCharts:
    def test_status_timeline_svg_present(self):
        sr = [_make_step_result(i) for i in range(1, 4)]
        html = make_dashboard_html("pipe", sr, 1.0, None, 5, 5)
        assert 'aria-label="Per-step status timeline"' in html

    def test_status_timeline_marks_cached_step(self):
        sr = [_make_step_result(1, cached=True)]
        html = make_dashboard_html("pipe", sr, 1.0, None, 5, 5)
        assert "replayed from cache" in html

    def test_duration_bars_svg_present(self):
        sr = [_make_step_result(i) for i in range(1, 4)]
        html = make_dashboard_html("pipe", sr, 1.0, None, 5, 5)
        assert 'aria-label="Per-step duration"' in html

    def test_slow_step_highlighted_in_warning_color(self):
        # Five fast steps plus one much slower one: the slow one should sit
        # at/above the 80th percentile and render in the warning color.
        sr = [_make_step_result(i, elapsed_sec=0.1) for i in range(1, 6)]
        sr.append(_make_step_result(6, elapsed_sec=50.0))
        html = make_dashboard_html("pipe", sr, 50.5, None, 5, 5)
        assert 'fill="var(--pc-warning)"' in html
        assert "(slow step)" in html

    def test_no_slow_highlight_when_all_steps_equal(self):
        sr = [_make_step_result(i, elapsed_sec=1.0) for i in range(1, 4)]
        html = make_dashboard_html("pipe", sr, 3.0, None, 5, 5)
        assert "(slow step)" not in html

    def test_site_flow_svg_present_with_legend(self):
        sr = [_make_step_result(i) for i in range(1, 4)]
        html = make_dashboard_html("pipe", sr, 1.0, None, 5, 5)
        assert 'aria-label="Site count in versus out per step"' in html
        assert "Sites in" in html
        assert "Sites out" in html

    def test_site_flow_end_labels_show_final_counts(self):
        sr = [_make_step_result(1), _make_step_result(2)]
        # override n_sites_out on the last step to a distinctive value
        sr[-1] = _make_step_result(2)
        sr[-1] = StepResult(
            step_idx=2,
            step_name="step_2",
            step_code="NR002",
            step_label="Step 2 Label",
            params={},
            elapsed_sec=0.2,
            n_sites_in=5,
            n_sites_out=3,
        )
        html = make_dashboard_html("pipe", sr, 1.0, None, 5, 3)
        assert ">3</text>" in html


# ─────────────────────────────────────────────────────────────────────────────
# Table view + step cards
# ─────────────────────────────────────────────────────────────────────────────


class TestTableAndCards:
    def test_table_view_present(self):
        sr = [_make_step_result(1)]
        html = make_dashboard_html("pipe", sr, 1.0, None, 5, 5)
        assert "<table>" in html
        assert "<th>Step</th>" in html

    def test_step_card_shows_error_box(self):
        sr = [_make_step_result(1, ok=False)]
        html = make_dashboard_html("pipe", sr, 0.5, None, 5, 0)
        assert 'class="error-box"' in html
        assert "test error" in html

    def test_step_card_omits_error_box_when_ok(self):
        sr = [_make_step_result(1, ok=True)]
        html = make_dashboard_html("pipe", sr, 0.5, None, 5, 5)
        assert 'class="error-box"' not in html

    def test_pipeline_yaml_block_included_when_given(self):
        sr = [_make_step_result(1)]
        html = make_dashboard_html(
            "pipe", sr, 0.5, None, 5, 5, pipeline_yaml="name: pipe\nsteps: []\n"
        )
        assert "Pipeline configuration" in html
        assert "name: pipe" in html

    def test_pipeline_yaml_block_omitted_when_empty(self):
        sr = [_make_step_result(1)]
        html = make_dashboard_html("pipe", sr, 0.5, None, 5, 5)
        assert "Pipeline configuration" not in html


# ─────────────────────────────────────────────────────────────────────────────
# Edge cases
# ─────────────────────────────────────────────────────────────────────────────


class TestEmptyStepResults:
    def test_empty_step_results_does_not_raise(self):
        html = make_dashboard_html("empty_pipe", [], 0.0, None, 0, 0)
        assert "<html" in html
        assert "empty_pipe" in html

    def test_empty_step_results_shows_zero_tiles(self):
        html = make_dashboard_html("empty_pipe", [], 0.0, None, 0, 0)
        assert "0/0 ok" in html

    def test_empty_step_results_charts_show_placeholder(self):
        html = make_dashboard_html("empty_pipe", [], 0.0, None, 0, 0)
        assert "No steps to display." in html


# ─────────────────────────────────────────────────────────────────────────────
# Real WILLY-data integration
# ─────────────────────────────────────────────────────────────────────────────

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
_DATA_DIR = _PROJECT_ROOT / "data" / "AMT" / "WILLY_DATA" / "L22PLT"
_OUTDIR = _PROJECT_ROOT / "pipe_results" / "dashboard_test"

pytestmark_willy = pytest.mark.skipif(
    not _DATA_DIR.exists(),
    reason=f"WILLY EDI data not found at {_DATA_DIR}",
)


@pytestmark_willy
class TestDashboardRealData:
    @pytest.fixture(scope="class")
    def dashboard_result(self):
        from pycsamt.emtools._core import ensure_sites
        from pycsamt.pipeline import Pipeline, Step, configure_pipe, reset_pipe

        configure_pipe(
            show_progress=False,
            plot_dpi=72,
            plot_fmt="png",
            report_formats=("html", "txt", "dashboard"),
        )

        if _OUTDIR.exists():
            shutil.rmtree(_OUTDIR)

        paths = sorted(_DATA_DIR.glob("*.edi"))[:5]
        sites = ensure_sites([str(p) for p in paths])

        pipe = Pipeline(
            [
                ("drop_dup", Step("FREQ002")),
                ("notch", Step("NR001", mains_hz=50)),
            ],
            name="willy_dashboard_demo",
        )
        result = pipe.run(
            sites,
            outdir=_OUTDIR,
            save_plots=True,
            save_edis=False,
            save_report=True,
        )
        yield result
        reset_pipe()

    def test_dashboard_html_written(self, dashboard_result):
        assert (_OUTDIR / "dashboard.html").is_file()

    def test_report_html_still_written(self, dashboard_result):
        # the fast tier keeps writing alongside the dashboard tier
        assert (_OUTDIR / "report.html").is_file()

    def test_dashboard_html_contains_pipeline_name(self, dashboard_result):
        html = (_OUTDIR / "dashboard.html").read_text(encoding="utf-8")
        assert "willy_dashboard_demo" in html

    def test_dashboard_html_contains_real_step_codes(self, dashboard_result):
        html = (_OUTDIR / "dashboard.html").read_text(encoding="utf-8")
        assert "FREQ002" in html
        assert "NR001" in html

    def test_dashboard_html_has_real_figures_count(self, dashboard_result):
        html = (_OUTDIR / "dashboard.html").read_text(encoding="utf-8")
        n_plots = len(dashboard_result.plots)
        assert f'<div class="value">{n_plots}</div>' in html
