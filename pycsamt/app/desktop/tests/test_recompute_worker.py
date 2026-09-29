# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for RecomputeWorker (pycsamt.app.desktop.workers.recompute_worker).

Strategy
--------
``EDIRecomputer`` (imported locally inside ``run()``) is monkeypatched at
its source module with a lightweight fake that invokes the worker's
``progress_callback`` like the real one does — this exercises the
worker's own signal-wiring and interruption-checking logic without
depending on real EDI I/O.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.workers.recompute_worker import RecomputeWorker


class _FakeResult:
    def __init__(self):
        self.n_ok = 2


def _fake_recomputer_cls(n_stations=2, raise_exc=None):
    class _FakeRecomputer:
        last_kwargs = None

        def __init__(self, **kw):
            _FakeRecomputer.last_kwargs = kw
            self._progress_callback = kw.get("progress_callback")

        def run(self, source):
            if raise_exc is not None:
                raise raise_exc
            if self._progress_callback:
                for i in range(1, n_stations + 1):
                    self._progress_callback(i, n_stations, f"STA{i}", "ok", "")
            return _FakeResult()

    return _FakeRecomputer


def _patch(monkeypatch, cls):
    import pycsamt.site.recompute as recompute_mod

    monkeypatch.setattr(recompute_mod, "EDIRecomputer", cls)


class TestSuccess:
    def test_run_success_emits_finished_and_progress(self, qapp, monkeypatch):
        _patch(monkeypatch, _fake_recomputer_cls(n_stations=3))
        w = RecomputeWorker(source=["fake"])
        results = []
        station_events = []
        w.finished.connect(results.append)
        w.station_done.connect(lambda *a: station_events.append(a))
        w.run()
        assert len(results) == 1
        assert w.result is results[0]
        assert len(station_events) == 3
        assert station_events[0] == (1, 3, "STA1", "ok", "")

    def test_constructor_forwards_kwargs(self, qapp, monkeypatch):
        fake_cls = _fake_recomputer_cls()
        _patch(monkeypatch, fake_cls)
        w = RecomputeWorker(
            source="src",
            output_root="/out",
            rotate_angle=15.0,
            fmin=1e-2,
            fmax=1e2,
            strict=True,
        )
        w.run()
        kw = fake_cls.last_kwargs
        assert kw["output_root"] == "/out"
        assert kw["rotate_angle"] == 15.0
        assert kw["fmin"] == pytest.approx(1e-2)
        assert kw["strict"] is True
        assert kw["progress"] is False
        assert kw["copy"] is True

    def test_result_attribute_defaults_none_before_run(self, qapp):
        w = RecomputeWorker(source=["fake"])
        assert w.result is None


class TestInterruption:
    def test_interruption_during_progress_callback_emits_error(self, qapp, monkeypatch):
        _patch(monkeypatch, _fake_recomputer_cls(n_stations=5))
        w = RecomputeWorker(source=["fake"])
        monkeypatch.setattr(w, "isInterruptionRequested", lambda: True)
        errors = []
        results = []
        w.error.connect(errors.append)
        w.finished.connect(results.append)
        w.run()
        assert results == []
        assert len(errors) == 1
        assert "cancelled" in errors[0].lower()


class TestExceptions:
    def test_recomputer_run_exception_emits_error(self, qapp, monkeypatch):
        _patch(
            monkeypatch,
            _fake_recomputer_cls(raise_exc=RuntimeError("recompute failed")),
        )
        w = RecomputeWorker(source=["fake"])
        errors = []
        results = []
        w.error.connect(errors.append)
        w.finished.connect(results.append)
        w.run()
        assert results == []
        assert errors == ["recompute failed"]
        assert w.result is None

    def test_recomputer_construction_exception_emits_error(self, qapp, monkeypatch):
        import pycsamt.site.recompute as recompute_mod

        def _boom(**kw):
            raise ValueError("bad config")

        monkeypatch.setattr(recompute_mod, "EDIRecomputer", _boom)
        w = RecomputeWorker(source=["fake"])
        errors = []
        w.error.connect(errors.append)
        w.run()
        assert errors == ["bad config"]
