# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for PhaseTensorStripGridDialog and _StripGridWorker
(pycsamt.app.desktop.tools.phase_tensor_strip_grid_tool).

Real data
---------
data/AMT/WILLY_DATA/L18PLT, L22PLT — real multi-line AMT EDIs. Station IDs
encode the survey line (``18-001A`` → line ``L18``, ``22-013VF`` → line
``L22``), which is exactly what ``_detect_lines`` needs to group on. A
5-station-per-line subset (10 EDIs total) is used for speed while still
giving ``_detect_lines`` two genuinely distinct real lines to group.
"""

from __future__ import annotations

import glob
from pathlib import Path
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.tools import phase_tensor_strip_grid_tool as tool_mod
from pycsamt.app.desktop.tools.phase_tensor_strip_grid_tool import (
    PhaseTensorStripGridDialog,
    _StripGridWorker,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY = _ROOT / "data" / "AMT" / "WILLY_DATA"
_L18 = _WILLY / "L18PLT"
_L22 = _WILLY / "L22PLT"

_HAS_WILLY = _L18.exists() and _L22.exists()


def _close():
    plt.close("all")


# ── Session-scoped fixtures ───────────────────────────────────────────────────


@pytest.fixture(scope="session")
def multi_line_sites():
    """10 real EDIs (5 from L18PLT + 5 from L22PLT) — two genuine lines."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY data not available")
    from pycsamt.emtools._core import ensure_sites

    l18 = sorted(glob.glob(str(_L18 / "*.edi")))[:5]
    l22 = sorted(glob.glob(str(_L22 / "*.edi")))[:5]
    return ensure_sites(l18 + l22)


@pytest.fixture(scope="session")
def single_line_sites():
    """5 real EDIs from L18PLT only — a single detected line."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY data not available")
    from pycsamt.emtools._core import ensure_sites

    l18 = sorted(glob.glob(str(_L18 / "*.edi")))[:5]
    return ensure_sites(l18)


# ── _detect_lines ─────────────────────────────────────────────────────────────


class TestDetectLines:
    def test_two_real_lines_detected(self, qapp, multi_line_sites):
        dlg = PhaseTensorStripGridDialog(sites=multi_line_sites)
        assert set(dlg._lines.keys()) == {"L18", "L22"}
        assert len(dlg._lines["L18"]) == 5
        assert len(dlg._lines["L22"]) == 5
        dlg.close()

    def test_lines_label_summarises_counts(self, qapp, multi_line_sites):
        dlg = PhaseTensorStripGridDialog(sites=multi_line_sites)
        text = dlg._lines_lbl.text()
        assert "L18" in text and "(5)" in text
        assert "L22" in text
        dlg.close()

    def test_single_line_detected(self, qapp, single_line_sites):
        dlg = PhaseTensorStripGridDialog(sites=single_line_sites)
        assert set(dlg._lines.keys()) == {"L18"}
        assert len(dlg._lines["L18"]) == 5
        dlg.close()

    def test_no_clear_grouping_falls_back_to_auto(self, qapp):
        """Station IDs with no separator all end up in their own group,
        which detect_lines_from_station_ids collapses to a single 'auto'
        bucket (see pycsamt.site.lines.detect_lines_from_station_ids)."""
        dlg = PhaseTensorStripGridDialog(sites=None)
        dlg._sites = [
            SimpleNamespace(station="AAA"),
            SimpleNamespace(station="BBB"),
            SimpleNamespace(station="CCC"),
        ]
        dlg._detect_lines()
        assert dlg._lines == {"auto": ["AAA", "BBB", "CCC"]}
        dlg.close()

    def test_no_sites_leaves_lines_empty(self, qapp):
        dlg = PhaseTensorStripGridDialog(sites=None)
        assert dlg._lines == {}
        assert dlg._lines_lbl.text() == "—"
        dlg.close()

    def test_unwrap_failure_per_station_is_swallowed(
        self, qapp, single_line_sites, monkeypatch
    ):
        """When _unwrap() raises for an item, _detect_lines must fall back
        to the raw (still-iterable) item rather than aborting the loop."""

        def _boom(_ed):
            raise RuntimeError("unwrap kaboom")

        monkeypatch.setattr("pycsamt.emtools._core._unwrap", _boom)
        dlg = PhaseTensorStripGridDialog(sites=single_line_sites)
        # The raw (non-unwrapped) items have no usable .station/.id, so
        # they all fall back to "?" and get grouped as a single line.
        assert dlg._lines
        dlg.close()

    def test_detect_lines_handles_exception(self, qapp, single_line_sites, monkeypatch):
        def _boom(_names, **_kw):
            raise RuntimeError("kaboom")

        monkeypatch.setattr("pycsamt.site.lines.detect_lines_from_station_ids", _boom)
        dlg = PhaseTensorStripGridDialog(sites=single_line_sites)
        assert dlg._lines == {}
        assert "failed" in dlg._lines_lbl.text().lower()
        assert "kaboom" in dlg._lines_lbl.text()
        dlg.close()


# ── _StripGridWorker ──────────────────────────────────────────────────────────


class TestStripGridWorker:
    def test_run_success_real_data(self, qapp, multi_line_sites):
        from pycsamt.emtools._core import _iter_items, _unwrap
        from pycsamt.site.lines import (
            detect_lines_from_station_ids,
            pick_representative_stations,
        )

        names = []
        for ed in _iter_items(multi_line_sites):
            try:
                ed = _unwrap(ed)
            except Exception:
                pass
            names.append(str(getattr(ed, "station", None) or "?"))
        lines = detect_lines_from_station_ids(names)
        profiles = {
            ln: pick_representative_stations(stns, 3) for ln, stns in lines.items()
        }

        done, error = [], []
        w = _StripGridWorker(multi_line_sites, profiles=profiles, c_by="skew")
        w.done.connect(done.append)
        w.error.connect(error.append)
        w.run()

        assert error == []
        assert len(done) == 1
        fig = done[0]
        assert hasattr(fig, "savefig")
        _close()

    def test_run_error_empty_profiles(self, qapp, multi_line_sites):
        done, error = [], []
        w = _StripGridWorker(multi_line_sites, profiles={}, c_by="skew")
        w.done.connect(done.append)
        w.error.connect(error.append)
        w.run()

        assert done == []
        assert len(error) == 1
        assert isinstance(error[0], str) and error[0]
        _close()

    def test_run_error_none_sites(self, qapp):
        done, error = [], []
        w = _StripGridWorker(None, profiles={"L18": ["18-001A"]}, c_by="skew")
        w.done.connect(done.append)
        w.error.connect(error.append)
        w.run()

        assert done == []
        assert len(error) == 1
        _close()


# ── Dialog construction ────────────────────────────────────────────────────────


class TestDialogConstruction:
    def test_construct_with_sites_enables_run(self, qapp, multi_line_sites):
        dlg = PhaseTensorStripGridDialog(sites=multi_line_sites)
        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == ""
        dlg.close()

    def test_construct_without_sites_disables_run(self, qapp):
        dlg = PhaseTensorStripGridDialog(sites=None)
        assert not dlg._run_btn.isEnabled()
        assert "No survey loaded." in dlg._status_lbl.text()
        dlg.close()

    def test_cby_combo_has_expected_items(self, qapp):
        dlg = PhaseTensorStripGridDialog(sites=None)
        items = [dlg._cby_combo.itemText(i) for i in range(dlg._cby_combo.count())]
        assert items == tool_mod._COLOR_BY
        dlg.close()

    def test_per_line_spin_default(self, qapp):
        dlg = PhaseTensorStripGridDialog(sites=None)
        assert dlg._per_line_spin.value() == 4
        dlg.close()

    def test_window_title_and_min_size(self, qapp):
        dlg = PhaseTensorStripGridDialog(sites=None)
        assert dlg.windowTitle() == "Phase Tensor Strip Grid"
        assert dlg.minimumWidth() == 1000
        assert dlg.minimumHeight() == 640
        dlg.close()


# ── _on_plot / _on_done / _on_error via fake worker ────────────────────────────


class _FakeSignal:
    def __init__(self):
        self.connected: list = []

    def connect(self, fn):
        self.connected.append(fn)

    def emit(self, *a, **k):
        for fn in list(self.connected):
            fn(*a, **k)


def _make_fake_worker(emit="done", payload=None):
    """Build a _StripGridWorker stand-in whose .start() synchronously
    emits either done(payload) or error(payload)."""
    captured: list = []

    class _Fake:
        def __init__(self, sites, profiles, c_by):
            self.sites = sites
            self.profiles = profiles
            self.c_by = c_by
            self.done = _FakeSignal()
            self.error = _FakeSignal()
            captured.append(self)

        def start(self):
            if emit == "done":
                self.done.emit(payload)
            else:
                self.error.emit(payload)

    _Fake.captured = captured
    return _Fake


class TestOnPlot:
    def test_no_lines_detected_shows_message(self, qapp, multi_line_sites):
        dlg = PhaseTensorStripGridDialog(sites=multi_line_sites)
        dlg._lines = {}
        dlg._on_plot()
        assert dlg._status_lbl.text() == "No survey lines detected."
        assert dlg._worker is None
        dlg.close()

    def test_builds_profiles_from_detected_lines_and_starts_worker(
        self, qapp, multi_line_sites, monkeypatch
    ):
        fake_fig = plt.figure()
        FakeWorker = _make_fake_worker(emit="done", payload=fake_fig)
        monkeypatch.setattr(tool_mod, "_StripGridWorker", FakeWorker)

        dlg = PhaseTensorStripGridDialog(sites=multi_line_sites)
        dlg._per_line_spin.setValue(3)
        dlg._cby_combo.setCurrentText("ellipt")
        dlg._on_plot()

        assert len(FakeWorker.captured) == 1
        w = FakeWorker.captured[0]
        assert w.sites is multi_line_sites
        assert w.c_by == "ellipt"
        assert set(w.profiles.keys()) == set(dlg._lines.keys())
        for _ln, stns in w.profiles.items():
            assert len(stns) <= 3
            assert stns == sorted(stns)

        # Fake worker's start() emitted done() synchronously -> _on_done ran.
        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Done."
        assert dlg._canvas.figure is fake_fig
        dlg.close()
        _close()

    def test_worker_error_path(self, qapp, multi_line_sites, monkeypatch):
        FakeWorker = _make_fake_worker(emit="error", payload="boom")
        monkeypatch.setattr(tool_mod, "_StripGridWorker", FakeWorker)

        dlg = PhaseTensorStripGridDialog(sites=multi_line_sites)
        dlg._on_plot()

        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Error: boom"
        dlg.close()

    def test_status_set_to_drawing_before_worker_starts(
        self, qapp, multi_line_sites, monkeypatch
    ):
        """Capture the status text set right before .start() is invoked by
        instrumenting a fake worker whose start() records the status seen
        so far (the real worker.start() is async, so this is the only
        directly-observable point for the 'Drawing…' message)."""
        seen = {}

        class _Fake:
            def __init__(self, sites, profiles, c_by):
                self.done = _FakeSignal()
                self.error = _FakeSignal()

            def start(self):
                seen["status"] = dlg._status_lbl.text()
                seen["run_enabled"] = dlg._run_btn.isEnabled()

        monkeypatch.setattr(tool_mod, "_StripGridWorker", _Fake)
        dlg = PhaseTensorStripGridDialog(sites=multi_line_sites)
        dlg._on_plot()

        assert seen["status"] == "Drawing ellipse-strip grid…"
        assert seen["run_enabled"] is False
        dlg.close()


class TestOnDoneOnError:
    def test_on_done_shows_figure_and_resets_status(self, qapp, multi_line_sites):
        dlg = PhaseTensorStripGridDialog(sites=multi_line_sites)
        dlg._run_btn.setEnabled(False)
        fig = plt.figure()
        dlg._on_done(fig)

        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Done."
        assert dlg._canvas.figure is fig
        dlg.close()
        _close()

    def test_on_error_sets_message_and_reenables_run(self, qapp, multi_line_sites):
        dlg = PhaseTensorStripGridDialog(sites=multi_line_sites)
        dlg._run_btn.setEnabled(False)
        dlg._on_error("something went wrong")

        assert dlg._run_btn.isEnabled()
        assert dlg._status_lbl.text() == "Error: something went wrong"
        dlg.close()
