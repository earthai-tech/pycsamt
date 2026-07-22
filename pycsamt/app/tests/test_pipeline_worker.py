# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for PipelineWorker (pycsamt.app.desktop.workers.pipeline_worker).

Strategy
--------
Driven against a lightweight fake controller (not the real
``PipelineController``) so each test can deterministically control which
steps succeed, fail, or are pre-marked SKIPPED, and exercise
``abort_on_error`` / interruption independent of which real pipeline
steps happen to be broken or slow at any given time.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers.pipeline_controller import StepStatus
from pycsamt.app.desktop.workers.pipeline_worker import PipelineWorker


class _FakeStep:
    def __init__(self, name, status=StepStatus.PENDING, abort_on_error=False):
        self.name = name
        self.status = status
        self.abort_on_error = abort_on_error
        self.elapsed_s = 0.0
        self.result_info = ""


class _FakeController:
    def __init__(self, steps, raises_at=None, results=None):
        self.steps = steps
        self._raises_at = raises_at or {}
        self._results = results or {}
        self.calls = []

    def execute_step(self, step_id, log_cb=None):
        self.calls.append(step_id)
        if log_cb:
            log_cb(f"executing step {step_id}")
        if step_id in self._raises_at:
            raise self._raises_at[step_id]
        self.steps[step_id].elapsed_s = 0.05
        self.steps[step_id].result_info = f"step {step_id} ok"
        return self._results.get(step_id, f"result-{step_id}")


def _run_worker(ctrl, step_ids):
    w = PipelineWorker(ctrl, step_ids)
    events = {
        "started": [],
        "done": [],
        "error": [],
        "skipped": [],
        "log": [],
        "progress": [],
        "all_done": [],
    }
    w.step_started.connect(events["started"].append)
    w.step_done.connect(lambda i, r: events["done"].append((i, r)))
    w.step_error.connect(lambda i, m: events["error"].append((i, m)))
    w.step_skipped.connect(events["skipped"].append)
    w.log_line.connect(events["log"].append)
    w.progress.connect(events["progress"].append)
    w.all_done.connect(lambda: events["all_done"].append(1))
    w.run()
    return events


class TestSuccessSequence:
    def test_all_steps_succeed(self, qapp):
        steps = [_FakeStep(f"Step{i}") for i in range(3)]
        ctrl = _FakeController(steps)
        events = _run_worker(ctrl, [0, 1, 2])

        assert events["started"] == [0, 1, 2]
        assert [i for i, _ in events["done"]] == [0, 1, 2]
        assert events["error"] == []
        assert events["progress"][-1] == 100
        assert events["all_done"] == [1]
        assert any("Step0" in ln for ln in events["log"])

    def test_progress_scales_with_step_count(self, qapp):
        steps = [_FakeStep(f"S{i}") for i in range(4)]
        ctrl = _FakeController(steps)
        events = _run_worker(ctrl, [0, 1, 2, 3])
        assert events["progress"] == [25, 50, 75, 100]


class TestSkippedSteps:
    def test_pre_skipped_step_emits_skipped_not_started(self, qapp):
        steps = [
            _FakeStep("A", status=StepStatus.SKIPPED),
            _FakeStep("B"),
        ]
        ctrl = _FakeController(steps)
        events = _run_worker(ctrl, [0, 1])

        assert events["skipped"] == [0]
        assert events["started"] == [1]
        assert 0 not in ctrl.calls
        assert any("skipped" in ln.lower() for ln in events["log"])
        assert events["progress"][-1] == 100


class TestErrorHandling:
    def test_step_error_without_abort_continues(self, qapp):
        steps = [
            _FakeStep("A"),
            _FakeStep("B"),
        ]
        ctrl = _FakeController(steps, raises_at={0: RuntimeError("boom")})
        events = _run_worker(ctrl, [0, 1])

        assert events["error"] == [(0, "boom")]
        assert [i for i, _ in events["done"]] == [1]
        assert events["all_done"] == [1]
        assert any("FAILED" in ln for ln in events["log"])

    def test_step_error_with_abort_stops_pipeline(self, qapp):
        steps = [
            _FakeStep("A", abort_on_error=True),
            _FakeStep("B"),
        ]
        ctrl = _FakeController(steps, raises_at={0: RuntimeError("fatal")})
        events = _run_worker(ctrl, [0, 1])

        assert events["error"] == [(0, "fatal")]
        assert events["done"] == []
        assert 1 not in ctrl.calls
        assert any("Aborting pipeline" in ln for ln in events["log"])
        assert events["all_done"] == [1]  # still emitted after the break


class TestInterruption:
    def test_interruption_requested_stops_before_first_step(
        self, qapp, monkeypatch
    ):
        steps = [_FakeStep("A"), _FakeStep("B")]
        ctrl = _FakeController(steps)
        w = PipelineWorker(ctrl, [0, 1])
        # QThread.requestInterruption() only takes effect on an actually
        # running thread; calling .run() synchronously here (as this whole
        # test suite does for determinism) never starts one, so the flag
        # never flips. Force it directly instead.
        monkeypatch.setattr(w, "isInterruptionRequested", lambda: True)
        events = {"started": [], "all_done": [], "log": []}
        w.step_started.connect(events["started"].append)
        w.all_done.connect(lambda: events["all_done"].append(1))
        w.log_line.connect(events["log"].append)
        w.run()

        assert events["started"] == []
        assert ctrl.calls == []
        assert events["all_done"] == [1]
        assert any("interrupted" in ln.lower() for ln in events["log"])


class TestLogFormat:
    def test_log_lines_are_timestamped(self, qapp):
        steps = [_FakeStep("A")]
        ctrl = _FakeController(steps)
        events = _run_worker(ctrl, [0])
        # "[HH:MM:SS]  message"
        assert all(ln.startswith("[") and "]" in ln for ln in events["log"])
