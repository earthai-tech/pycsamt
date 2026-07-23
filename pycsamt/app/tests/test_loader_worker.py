# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for LoaderWorker (pycsamt.app.desktop.workers.loader_worker).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — a handful of real EDI files loaded through
the real ``DataController.load()`` for the success path, since EDI
parsing is fast and deterministic (no mocking needed). The failure path
is exercised with a nonexistent path.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers.data_controller import DataController
from pycsamt.app.desktop.workers.loader_worker import LoaderWorker

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


@pytest.fixture
def willy_paths():
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    return [str(p) for p in sorted(_WILLY_L18.glob("*.edi"))[:3]]


class TestLoaderWorker:
    def test_construction_builds_data_controller(self, qapp):
        w = LoaderWorker(paths=["/some/path.edi"])
        assert isinstance(w.data_controller, DataController)

    def test_load_success_emits_finished(self, qapp, willy_paths):
        w = LoaderWorker(paths=willy_paths)
        results = []
        errors = []
        progresses = []
        w.finished.connect(results.append)
        w.error.connect(errors.append)
        w.progress.connect(progresses.append)
        w.run()
        assert errors == []
        assert len(results) == 1
        assert len(progresses) > 0  # progress_callback wired through

    def test_load_failure_emits_error(self, qapp, monkeypatch):
        w = LoaderWorker(paths=["/definitely/not/a/real/path.edi"])

        def _boom(paths):
            raise RuntimeError("cannot load")

        monkeypatch.setattr(w._ctrl, "load", _boom)
        results = []
        errors = []
        w.finished.connect(results.append)
        w.error.connect(errors.append)
        w.run()
        assert results == []
        assert errors == ["cannot load"]

    def test_data_controller_property_returns_same_instance_after_load(
        self, qapp, willy_paths
    ):
        w = LoaderWorker(paths=willy_paths)
        ctrl_before = w.data_controller
        w.run()
        assert w.data_controller is ctrl_before
