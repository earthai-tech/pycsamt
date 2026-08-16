"""Tests for chunk-4 observability: on_step callback, progress_style="rich",
and PipelineResult._repr_html_.

Public entry point: :mod:`pycsamt.pipeline` (``Pipeline.run(on_step=...)``,
``configure_pipe(progress_style="rich")``, ``PipelineResult._repr_html_``)
Implementation: :mod:`pycsamt.pipeline._pipeline`, :mod:`pycsamt.pipeline._richview`,
:mod:`pycsamt.pipeline._report`

Test doubles mirror the ones in ``test_pipeline.py`` (``_IdentityStep``,
``_MutatingStep``, ``_CountingSites``, ``_raising``) but are redefined
locally here, the same self-contained-per-file convention
``test_cache.py`` already established in chunk 3.
"""

from __future__ import annotations

import warnings
from typing import Any

import pytest

from pycsamt.api.pipe import configure_pipe, reset_pipe
from pycsamt.pipeline import Pipeline, Step, StepResult, lookup_step


class _CountingSites:
    def __init__(self, n: int = 3, count: int = 0):
        self._n = n
        self.count = count

    def __len__(self) -> int:
        return self._n


class _MutatingStep(Step):
    def transform(self, sites):
        return _CountingSites(n=sites._n, count=sites.count + 1)

    def generate_qc_plots(self, sites):
        return []


def _mutating(code: str) -> _MutatingStep:
    step = _MutatingStep.__new__(_MutatingStep)
    step.spec = lookup_step(code)
    step.params = dict(step.spec.defaults)
    return step


class _RaisingStep(Step):
    def transform(self, sites):
        raise RuntimeError("deliberate step failure")

    def generate_qc_plots(self, sites):
        return []


def _raising(code: str) -> _RaisingStep:
    step = _RaisingStep.__new__(_RaisingStep)
    step.spec = lookup_step(code)
    step.params = dict(step.spec.defaults)
    return step


@pytest.fixture(autouse=True)
def _reset_pipe_cfg():
    reset_pipe()
    yield
    reset_pipe()


@pytest.fixture()
def sites() -> _CountingSites:
    return _CountingSites(n=5, count=0)


@pytest.fixture()
def simple_pipe() -> Pipeline:
    return Pipeline(
        [("a", _mutating("NR001")), ("b", _mutating("FREQ001"))],
        name="obs_test",
    )


# ─────────────────────────────────────────────────────────────────────────────
# on_step
# ─────────────────────────────────────────────────────────────────────────────


class TestOnStepCallback:
    def test_fires_once_per_step_in_order(self, simple_pipe, sites):
        seen: list[str] = []
        simple_pipe.run(
            sites,
            outdir=None,
            save_report=False,
            save_plots=False,
            on_step=lambda sr: seen.append(sr.step_code),
        )
        assert seen == ["NR001", "FREQ001"]

    def test_receives_real_stepresult(self, simple_pipe, sites):
        received: list[StepResult] = []
        simple_pipe.run(
            sites,
            outdir=None,
            save_report=False,
            save_plots=False,
            on_step=received.append,
        )
        assert all(isinstance(sr, StepResult) for sr in received)
        assert received[0].step_idx == 1
        assert received[1].step_idx == 2
        # _CountingSites.__len__ always reports _n (unchanged by mutation);
        # real proof of chaining is in TestPipelineRun.test_run_chaining_increments_count
        assert received[0].n_sites_in == sites._n
        assert received[1].n_sites_in == sites._n

    def test_none_is_default_and_is_a_no_op(self, simple_pipe, sites):
        result = simple_pipe.run(sites, outdir=None, save_report=False, save_plots=False)
        assert result.ok  # no crash, unchanged behavior

    def test_callback_exception_does_not_abort_run(self, simple_pipe, sites):
        def bad(sr):
            raise RuntimeError("callback boom")

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = simple_pipe.run(
                sites, outdir=None, save_report=False, save_plots=False, on_step=bad
            )
        assert result.ok
        assert len(result.step_results) == 2
        assert any("callback boom" in str(w.message) for w in caught)

    def test_does_not_fire_for_a_step_that_raises_under_raise_policy(self, sites):
        configure_pipe(on_step_error="raise")
        pipe = Pipeline([("bad", _raising("NR001")), ("ok", _mutating("FREQ001"))])
        seen: list[str] = []
        with pytest.raises(RuntimeError, match="deliberate step failure"):
            pipe.run(
                sites,
                outdir=None,
                save_report=False,
                save_plots=False,
                on_step=lambda sr: seen.append(sr.step_code),
            )
        assert seen == []  # the failing step never got a StepResult

    def test_fires_for_warn_policy_error_step_too(self, sites):
        configure_pipe(on_step_error="warn")
        pipe = Pipeline([("bad", _raising("NR001")), ("ok", _mutating("FREQ001"))])
        seen: list[tuple[str, bool]] = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pipe.run(
                sites,
                outdir=None,
                save_report=False,
                save_plots=False,
                on_step=lambda sr: seen.append((sr.step_code, sr.ok)),
            )
        assert seen == [("NR001", False), ("FREQ001", True)]

    def test_fires_on_cache_hit_with_cached_true(self, simple_pipe, sites, tmp_path):
        simple_pipe.run(sites, outdir=None, save_report=False, save_plots=False, cache=tmp_path)
        seen: list[bool] = []
        simple_pipe.run(
            sites,
            outdir=None,
            save_report=False,
            save_plots=False,
            cache=tmp_path,
            on_step=lambda sr: seen.append(sr.cached),
        )
        assert seen == [True, True]


# ─────────────────────────────────────────────────────────────────────────────
# progress_style="rich"
# ─────────────────────────────────────────────────────────────────────────────


class TestRichProgressStyle:
    def test_run_completes_successfully(self, simple_pipe, sites):
        configure_pipe(progress_style="rich", show_progress=True)
        result = simple_pipe.run(sites, outdir=None, save_report=False, save_plots=False)
        assert result.ok
        assert len(result.step_results) == 2

    def test_rich_style_does_not_use_log_style_markers(self, simple_pipe, sites, capsys):
        configure_pipe(progress_style="rich", show_progress=True)
        simple_pipe.run(sites, outdir=None, save_report=False, save_plots=False)
        out = capsys.readouterr().out
        # "log" style's distinctive per-step marker must not appear when
        # progress_style="rich" is selected instead.
        assert "[1/2] a [NR001]" not in out
        # the rich table's own content should appear somewhere in plain text
        assert "NR001" in out

    def test_silent_when_show_progress_false(self, simple_pipe, sites, capsys):
        configure_pipe(progress_style="rich", show_progress=False)
        simple_pipe.run(sites, outdir=None, save_report=False, save_plots=False)
        out = capsys.readouterr().out
        assert out == ""

    def test_works_with_cache_and_on_step_simultaneously(self, simple_pipe, sites, tmp_path):
        configure_pipe(progress_style="rich", show_progress=True)
        seen = []
        result = simple_pipe.run(
            sites,
            outdir=None,
            save_report=False,
            save_plots=False,
            cache=tmp_path,
            on_step=seen.append,
        )
        assert result.ok
        assert len(seen) == 2


# ─────────────────────────────────────────────────────────────────────────────
# PipelineResult._repr_html_
# ─────────────────────────────────────────────────────────────────────────────


class TestReprHtml:
    def test_returns_html_string(self, simple_pipe, sites):
        result = simple_pipe.run(sites, outdir=None, save_report=False, save_plots=False)
        html = result._repr_html_()
        assert isinstance(html, str)
        assert "<table>" in html
        assert "</table>" in html

    def test_contains_pipeline_name_and_step_codes(self, simple_pipe, sites):
        result = simple_pipe.run(sites, outdir=None, save_report=False, save_plots=False)
        html = result._repr_html_()
        assert "obs_test" in html
        assert "NR001" in html
        assert "FREQ001" in html

    def test_shows_cached_marker(self, simple_pipe, sites, tmp_path):
        simple_pipe.run(sites, outdir=None, save_report=False, save_plots=False, cache=tmp_path)
        result = simple_pipe.run(
            sites, outdir=None, save_report=False, save_plots=False, cache=tmp_path
        )
        html = result._repr_html_()
        assert "cached" in html.lower()

    def test_handles_empty_step_results(self):
        from pycsamt.pipeline import PipelineResult

        result = PipelineResult(
            sites_in=None, sites_out=None, step_results=[], outdir=None, elapsed_sec=0.0
        )
        html = result._repr_html_()
        assert isinstance(html, str)
        assert "<table>" in html

    def test_scoped_css_does_not_leak_global_selectors(self, simple_pipe, sites):
        """A real bug class this must avoid: injecting unscoped `table {...}`
        CSS into a notebook cell would restyle every other output on the
        page. Every rule must be scoped under .pycsamt-result."""
        result = simple_pipe.run(sites, outdir=None, save_report=False, save_plots=False)
        html = result._repr_html_()
        style_block = html.split("<style>")[1].split("</style>")[0]
        for line in style_block.splitlines():
            line = line.strip()
            if "{" not in line:
                continue
            selector = line.split("{")[0].strip()
            assert selector.startswith(".pycsamt-result"), (
                f"unscoped CSS selector would leak into the notebook page: {selector!r}"
            )
