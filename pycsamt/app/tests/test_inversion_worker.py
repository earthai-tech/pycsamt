# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for InversionWorker (pycsamt.app.desktop.workers.inversion_worker).

Strategy
--------
A real Occam2D run needs a compiled binary and real input files, so
``OccamRunner`` / ``InversionResult`` (imported locally inside
``_run_inversion``) are monkeypatched at their source module
(``pycsamt.models.occam2d``) with lightweight fakes that preserve the
real constructor/``run()`` contract, including the optional
``iter_callback`` attribute the worker probes for via ``hasattr``.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.workers.inversion_worker import InversionWorker


class _FakeRunnerNoCallback:
    """Fake OccamRunner that does NOT support per-iteration callbacks."""

    def __init__(self, workdir, binary_path):
        self.workdir = workdir
        self.binary_path = binary_path

    def run(self, max_iter, target_misfit):
        return 0


class _FakeRunnerWithCallback:
    """Fake OccamRunner that invokes iter_callback like the real one."""

    def __init__(self, workdir, binary_path):
        self.workdir = workdir
        self.binary_path = binary_path
        self.iter_callback = None

    def run(self, max_iter, target_misfit):
        for i in range(1, 4):
            self.iter_callback(i, 5.0 / i)
        return 0


class _FakeInversionResult:
    def __init__(self, workdir):
        self.workdir = workdir
        self.summary = f"fake summary for {workdir}"


def _patch_occam(monkeypatch, runner_cls, result_cls=_FakeInversionResult):
    import pycsamt.models.occam2d as occam_mod

    monkeypatch.setattr(occam_mod, "OccamRunner", runner_cls)
    monkeypatch.setattr(occam_mod, "InversionResult", result_cls)


@pytest.fixture
def workdir(tmp_path):
    d = tmp_path / "occam_run"
    d.mkdir()
    return d


class TestGuards:
    def test_missing_workdir_reports_error(self, qapp, tmp_path):
        w = InversionWorker(workdir=str(tmp_path / "does_not_exist"))
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert len(errors) == 1
        assert "not found" in errors[0]

    def test_missing_binary_reports_error(self, qapp, workdir):
        w = InversionWorker(
            workdir=str(workdir), binary_path=str(workdir / "no_such_binary")
        )
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert len(errors) == 1
        assert "binary not found" in errors[0]


class TestSuccess:
    def test_run_without_iter_callback(self, qapp, workdir, monkeypatch):
        _patch_occam(monkeypatch, _FakeRunnerNoCallback)
        w = InversionWorker(workdir=str(workdir), max_iter=10, target_misfit=2.0)
        results = []
        lines = []
        progresses = []
        w.finished.connect(results.append)
        w.stdout_line.connect(lines.append)
        w.progress.connect(progresses.append)
        w.run()
        assert len(results) == 1
        assert progresses[-1] == 100
        assert any("output after completion" in ln for ln in lines)
        assert any("fake summary" in ln for ln in lines)

    def test_run_with_iter_callback_emits_per_iteration(
        self, qapp, workdir, monkeypatch
    ):
        _patch_occam(monkeypatch, _FakeRunnerWithCallback)
        w = InversionWorker(workdir=str(workdir), max_iter=3, target_misfit=1.0)
        results = []
        lines = []
        w.finished.connect(results.append)
        w.stdout_line.connect(lines.append)
        w.run()
        assert len(results) == 1
        assert any("iter   1" in ln for ln in lines)
        assert any("iter   3" in ln for ln in lines)

    def test_binary_provided_and_valid(self, qapp, workdir, monkeypatch):
        _patch_occam(monkeypatch, _FakeRunnerNoCallback)
        binary = workdir / "occam2d_binary"
        binary.write_text("fake binary")
        w = InversionWorker(workdir=str(workdir), binary_path=str(binary))
        lines = []
        w.stdout_line.connect(lines.append)
        w.run()
        assert any(str(binary) in ln for ln in lines)


class TestCancellation:
    def test_cancel_sets_flag(self, qapp, workdir):
        w = InversionWorker(workdir=str(workdir))
        assert not w._cancelled
        w.cancel()
        assert w._cancelled

    def test_cancelled_without_callback_skips_result(
        self, qapp, workdir, monkeypatch
    ):
        _patch_occam(monkeypatch, _FakeRunnerNoCallback)
        w = InversionWorker(workdir=str(workdir))
        w.cancel()
        results = []
        lines = []
        w.finished.connect(results.append)
        w.stdout_line.connect(lines.append)
        w.run()
        assert results == []
        assert any("cancelled" in ln.lower() for ln in lines)

    def test_cancelled_mid_run_via_callback_reports_error(
        self, qapp, workdir, monkeypatch
    ):
        """
        With a runner that supports per-iteration callbacks, cancellation
        is detected *inside* the callback by raising InterruptedError --
        which propagates out of run() and is caught by the worker's own
        broad except, surfacing as an error rather than the clean
        post-run "cancelled" message (that branch is effectively
        unreachable once iteration has already started).
        """
        _patch_occam(monkeypatch, _FakeRunnerWithCallback)
        w = InversionWorker(workdir=str(workdir))
        w.cancel()
        results = []
        errors = []
        w.finished.connect(results.append)
        w.error.connect(errors.append)
        w.run()
        assert results == []
        assert len(errors) == 1
        assert "Cancelled by user" in errors[0]


class TestExceptions:
    def test_runner_construction_failure_reports_error(
        self, qapp, workdir, monkeypatch
    ):
        import pycsamt.models.occam2d as occam_mod

        def _boom(**kw):
            raise RuntimeError("cannot build runner")

        monkeypatch.setattr(occam_mod, "OccamRunner", _boom)
        w = InversionWorker(workdir=str(workdir))
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert errors == ["cannot build runner"]

    def test_runner_run_failure_reports_error(self, qapp, workdir, monkeypatch):
        class _BoomRunner:
            def __init__(self, workdir, binary_path):
                pass

            def run(self, max_iter, target_misfit):
                raise RuntimeError("solver crashed")

        _patch_occam(monkeypatch, _BoomRunner)
        w = InversionWorker(workdir=str(workdir))
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert errors == ["solver crashed"]
