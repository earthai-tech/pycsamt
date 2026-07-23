# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for FormatConverterDialog and _ConvertWorker
(pycsamt.app.desktop.tools.converter_tool).

Strategy
--------
* ``_ConvertWorker`` is a real ``QThread`` subclass. Its ``.run()`` method
  is called directly (synchronously, no real threading) against small
  real ``Sites`` collections built from a handful of real EDI files, so
  the actual metadata-extraction and file-writing logic is exercised —
  mirroring the ``InversionWorker.run()`` idiom used in
  ``test_inversion_dlg.py``.
* ``FormatConverterDialog``'s own orchestration is tested by monkeypatching
  ``converter_tool._ConvertWorker`` with a lightweight fake class whose
  ``.start()`` synchronously emits through plain-callable "signals"
  (``_FakeSignal``) — the same idiom used in ``test_main_window_actions.py``
  and ``test_recompute_dlg.py`` — so no real QThread is ever spun up.

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — first 3 EDI files (small subset, for speed).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.tools import converter_tool
from pycsamt.app.desktop.tools.converter_tool import (
    FormatConverterDialog,
    _ConvertWorker,
)

# ── Real data ────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # repo root
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


@pytest.fixture(scope="module")
def small_sites():
    """3-station Sites collection loaded once for this module."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    files = sorted(_WILLY_L18.glob("*.edi"))[:3]
    return ensure_sites([str(f) for f in files])


# ── Fake worker / signal idiom (mirrors test_recompute_dlg.py) ───────────────


class _FakeSignal:
    def __init__(self):
        self._fns = []

    def connect(self, fn):
        self._fns.append(fn)

    def emit(self, *a):
        for fn in list(self._fns):
            fn(*a)


def _fake_worker_cls(
    *, progress_events=(), done_msg=None, error_msg=None, on_start=None
):
    """Build a lightweight _ConvertWorker stand-in.

    ``.start()`` synchronously replays ``progress_events`` then emits
    either ``done`` or ``error`` — no real QThread involved.
    """

    class _FakeWorker:
        captured = []

        def __init__(self, sites, out_dir, fmt):
            self.sites = sites
            self.out_dir = out_dir
            self.fmt = fmt
            self.progress = _FakeSignal()
            self.done = _FakeSignal()
            self.error = _FakeSignal()
            _FakeWorker.captured.append(self)

        def start(self):
            if on_start is not None:
                on_start(self)
            for ev in progress_events:
                self.progress.emit(*ev)
            if error_msg is not None:
                self.error.emit(error_msg)
            elif done_msg is not None:
                self.done.emit(done_msg)

    return _FakeWorker


# ── _ConvertWorker — CSV format ───────────────────────────────────────────────


class TestConvertWorkerCSV:
    def test_writes_csv_file(self, tmp_path, small_sites):
        worker = _ConvertWorker(
            small_sites, tmp_path, "CSV (station metadata)"
        )
        done, errors = [], []
        worker.done.connect(done.append)
        worker.error.connect(errors.append)
        worker.run()

        assert errors == []
        assert len(done) == 1
        out_path = tmp_path / "survey_stations.csv"
        assert out_path.exists()
        text = out_path.read_text(encoding="utf-8")
        assert "station" in text.splitlines()[0]

    def test_progress_emitted_per_station(self, tmp_path, small_sites):
        worker = _ConvertWorker(
            small_sites, tmp_path, "CSV (station metadata)"
        )
        calls = []
        worker.progress.connect(lambda cur, total, name: calls.append(
            (cur, total, name)
        ))
        worker.run()

        assert len(calls) == 3
        assert calls[-1][0] == 3  # last "cur" == total
        assert all(c[1] == 3 for c in calls)  # total is constant

    def test_done_message_mentions_count_and_dir(self, tmp_path, small_sites):
        worker = _ConvertWorker(
            small_sites, tmp_path, "CSV (station metadata)"
        )
        done = []
        worker.done.connect(done.append)
        worker.run()

        assert "3" in done[0]
        assert str(tmp_path) in done[0]

    def test_empty_sites_raises_write_error(self, tmp_path):
        """rows is empty -> rows[0] in the CSV branch raises IndexError,
        which is caught and re-emitted as an `error` signal."""
        worker = _ConvertWorker([], tmp_path, "CSV (station metadata)")
        done, errors = [], []
        worker.done.connect(done.append)
        worker.error.connect(errors.append)
        worker.run()

        assert done == []
        assert len(errors) == 1
        assert "Write error" in errors[0]


# ── _ConvertWorker — JSON format ──────────────────────────────────────────────


class TestConvertWorkerJSON:
    def test_writes_json_file(self, tmp_path, small_sites):
        worker = _ConvertWorker(
            small_sites, tmp_path, "JSON (station metadata)"
        )
        done, errors = [], []
        worker.done.connect(done.append)
        worker.error.connect(errors.append)
        worker.run()

        assert errors == []
        assert len(done) == 1
        out_path = tmp_path / "survey_stations.json"
        assert out_path.exists()
        rows = json.loads(out_path.read_text(encoding="utf-8"))
        assert isinstance(rows, list)
        assert len(rows) == 3
        for key in (
            "station",
            "lat",
            "lon",
            "n_freq",
            "t_min",
            "t_max",
            "has_z_err",
        ):
            assert key in rows[0]
        assert rows[0]["n_freq"] > 0

    def test_empty_sites_still_succeeds(self, tmp_path):
        """Unlike CSV, an empty rows list is valid JSON ('[]'), so no
        error is raised for the JSON branch."""
        worker = _ConvertWorker([], tmp_path, "JSON (station metadata)")
        done, errors = [], []
        worker.done.connect(done.append)
        worker.error.connect(errors.append)
        worker.run()

        assert errors == []
        assert len(done) == 1
        out_path = tmp_path / "survey_stations.json"
        assert json.loads(out_path.read_text(encoding="utf-8")) == []


# ── _ConvertWorker — EDI format ───────────────────────────────────────────────


class TestConvertWorkerEDI:
    def test_edi_format_writes_no_edi_files(self, tmp_path, small_sites):
        """BUG: the EDI re-export branch looks for a ``write_edi_file``
        method (``getattr(edi_obj, "write_edi_file", None)``), but no
        object in this codebase exposes that name — ``EDIFile`` exposes
        ``write`` instead (see pycsamt/site/base.py:1573). So choosing
        "EDI (re-export)" silently writes *zero* .edi files while still
        reporting success via the `done` signal. Documented here, not
        fixed (tests-only constraint)."""
        worker = _ConvertWorker(small_sites, tmp_path, "EDI (re-export)")
        done, errors = [], []
        worker.done.connect(done.append)
        worker.error.connect(errors.append)
        worker.run()

        assert errors == []
        assert len(done) == 1
        assert tmp_path.exists()
        assert list(tmp_path.iterdir()) == []

    def test_edi_format_creates_out_dir(self, tmp_path, small_sites):
        out_dir = tmp_path / "nested" / "export"
        worker = _ConvertWorker(small_sites, out_dir, "EDI (re-export)")
        worker.run()
        assert out_dir.exists()


# ── _ConvertWorker — metadata edge cases (synthetic fake EDI-like items) ─────


class _FakeZNoFreq:
    freq = None
    z_err = None


class _FakeZAllNonPositiveFreq:
    freq = [0.0, -1.0]
    z_err = None


class _FakeEdInvalidLatLon:
    station = "BAD1"
    lat = "not-a-number"
    lon = "not-a-number"
    Z = None


class _FakeEdNoZFreq:
    station = "NOZ1"
    lat = 1.0
    lon = 2.0
    Z = _FakeZNoFreq()


class _FakeEdNonPositiveFreq:
    station = "NPF1"
    lat = 1.0
    lon = 2.0
    Z = _FakeZAllNonPositiveFreq()


class _FakeEdiObjWithWriter:
    def __init__(self):
        self.calls = []

    def write_edi_file(self, path):
        self.calls.append(path)


class _FakeEdWithWriter:
    station = "W1"
    lat = 0.0
    lon = 0.0
    Z = None

    def __init__(self):
        self.edi = _FakeEdiObjWithWriter()


class TestConvertWorkerMetadataEdgeCases:
    def test_invalid_lat_lon_written_as_nan(self, tmp_path):
        worker = _ConvertWorker(
            [_FakeEdInvalidLatLon()], tmp_path, "JSON (station metadata)"
        )
        done = []
        worker.done.connect(done.append)
        worker.run()

        assert len(done) == 1
        rows = json.loads(
            (tmp_path / "survey_stations.json").read_text(encoding="utf-8")
        )
        assert rows[0]["lat"] != rows[0]["lat"]  # NaN
        assert rows[0]["lon"] != rows[0]["lon"]  # NaN

    def test_z_object_with_no_freq_yields_zero_n_freq(self, tmp_path):
        worker = _ConvertWorker(
            [_FakeEdNoZFreq()], tmp_path, "JSON (station metadata)"
        )
        worker.run()
        rows = json.loads(
            (tmp_path / "survey_stations.json").read_text(encoding="utf-8")
        )
        assert rows[0]["n_freq"] == 0
        assert rows[0]["t_min"] is None
        assert rows[0]["t_max"] is None

    def test_z_object_with_only_non_positive_freq_yields_zero_n_freq(
        self, tmp_path
    ):
        worker = _ConvertWorker(
            [_FakeEdNonPositiveFreq()], tmp_path, "JSON (station metadata)"
        )
        worker.run()
        rows = json.loads(
            (tmp_path / "survey_stations.json").read_text(encoding="utf-8")
        )
        assert rows[0]["n_freq"] == 0

    def test_write_edi_file_is_invoked_when_present(self, tmp_path):
        """Counterpart to the EDI-bug test: when an item *does* expose a
        ``write_edi_file`` method, the worker does call it."""
        fake = _FakeEdWithWriter()
        worker = _ConvertWorker([fake], tmp_path, "EDI (re-export)")
        done = []
        worker.done.connect(done.append)
        worker.run()

        assert len(done) == 1
        assert fake.edi.calls, "write_edi_file should have been invoked"

    def test_write_edi_file_exception_is_swallowed(self, tmp_path):
        """If the item's `write_edi_file` itself raises, the worker must
        catch it and keep going (best-effort export)."""

        class _RaisingEdiObj:
            def write_edi_file(self, path):
                raise OSError("disk full")

        class _FakeEdRaisingWriter:
            station = "R1"
            lat = 0.0
            lon = 0.0
            Z = None

            def __init__(self):
                self.edi = _RaisingEdiObj()

        worker = _ConvertWorker(
            [_FakeEdRaisingWriter()], tmp_path, "EDI (re-export)"
        )
        done, errors = [], []
        worker.done.connect(done.append)
        worker.error.connect(errors.append)
        worker.run()

        assert errors == []
        assert len(done) == 1

    def test_unwrap_exception_is_swallowed(
        self, monkeypatch, tmp_path, small_sites
    ):
        """When `_unwrap` itself raises, the worker must catch it and
        keep processing the original (non-unwrapped) item."""
        import pycsamt.emtools._core as core

        def _raise(_ed):
            raise RuntimeError("boom")

        monkeypatch.setattr(core, "_unwrap", _raise)

        worker = _ConvertWorker(
            small_sites, tmp_path, "JSON (station metadata)"
        )
        done, errors = [], []
        worker.done.connect(done.append)
        worker.error.connect(errors.append)
        worker.run()

        assert errors == []
        assert len(done) == 1


# ── _ConvertWorker — broken _core import fallback path ───────────────────────


class TestConvertWorkerImportFallback:
    def test_core_import_failure_falls_back_to_plain_iteration(
        self, monkeypatch, tmp_path
    ):
        """When `from pycsamt.emtools._core import _iter_items, _unwrap`
        itself fails (e.g. a broken installation), the worker falls back
        to plain ``list(self._sites)`` iteration instead of crashing."""
        import sys
        import types

        dummy = types.ModuleType("pycsamt.emtools._core")
        monkeypatch.setitem(sys.modules, "pycsamt.emtools._core", dummy)

        worker = _ConvertWorker(
            ["a", "b"], tmp_path, "JSON (station metadata)"
        )
        done, errors = [], []
        worker.done.connect(done.append)
        worker.error.connect(errors.append)
        worker.run()

        assert errors == []
        assert len(done) == 1
        rows = json.loads(
            (tmp_path / "survey_stations.json").read_text(encoding="utf-8")
        )
        assert len(rows) == 2

    def test_core_import_failure_and_non_iterable_sites_emits_error(
        self, monkeypatch, tmp_path
    ):
        """Both the primary and fallback iteration paths fail -> the
        worker must emit `error`, not raise."""
        import sys
        import types

        dummy = types.ModuleType("pycsamt.emtools._core")
        monkeypatch.setitem(sys.modules, "pycsamt.emtools._core", dummy)

        class _BadIterable:
            def __iter__(self):
                raise RuntimeError("cannot iterate")

        worker = _ConvertWorker(
            _BadIterable(), tmp_path, "JSON (station metadata)"
        )
        done, errors = [], []
        worker.done.connect(done.append)
        worker.error.connect(errors.append)
        worker.run()

        assert done == []
        assert len(errors) == 1
        assert "Cannot iterate survey" in errors[0]


# ── _ConvertWorker — error path (invalid out_dir) ─────────────────────────────


class TestConvertWorkerMkdirFailure:
    def test_mkdir_failure_is_not_caught(self, tmp_path):
        """BUG: ``self._out_dir.mkdir(parents=True, exist_ok=True)`` is
        called outside any try/except. When ``out_dir``'s parent is
        actually a *file* (not a directory), mkdir raises an OSError
        uncaught (``NotADirectoryError`` on POSIX, ``FileExistsError``
        on Windows) — neither `done` nor `error` is ever emitted. In
        real QThread usage (``.start()``) this means the exception is
        silently swallowed by Qt and the dialog's Convert button stays
        disabled forever, with no feedback to the user. Documented
        here, not fixed."""
        blocker = tmp_path / "blocker.txt"
        blocker.write_text("x", encoding="utf-8")
        out_dir = blocker / "sub"

        worker = _ConvertWorker([], out_dir, "CSV (station metadata)")
        with pytest.raises(OSError):
            worker.run()


# ── FormatConverterDialog — construction ──────────────────────────────────────


class TestDialogConstruction:
    def test_creates_with_no_sites(self, qapp):
        dlg = FormatConverterDialog(None)
        assert dlg.windowTitle() == "Format Converter"
        assert dlg._src_lbl.text() == "No survey loaded"
        dlg.close()

    def test_creates_with_sites_shows_count(self, qapp, small_sites):
        dlg = FormatConverterDialog(small_sites)
        assert "3" in dlg._src_lbl.text()
        assert "Loaded survey" in dlg._src_lbl.text()
        dlg.close()

    def test_format_combo_has_three_options(self, qapp):
        dlg = FormatConverterDialog(None)
        items = [
            dlg._fmt_combo.itemText(i)
            for i in range(dlg._fmt_combo.count())
        ]
        assert items == [
            "CSV (station metadata)",
            "JSON (station metadata)",
            "EDI (re-export)",
        ]
        dlg.close()

    def test_default_out_dir_is_home_export(self, qapp):
        dlg = FormatConverterDialog(None)
        assert dlg._dir_edit.text() == str(Path.home() / "pycsamt_export")
        dlg.close()

    def test_progress_starts_at_zero(self, qapp):
        dlg = FormatConverterDialog(None)
        assert dlg._progress.value() == 0
        dlg.close()

    def test_build_ui_handles_core_import_failure(
        self, qapp, monkeypatch, small_sites
    ):
        """When `_iter_items` can't be imported while building the
        source label, the dialog must fall back to the plain
        "Loaded survey" label instead of raising."""
        import sys
        import types

        dummy = types.ModuleType("pycsamt.emtools._core")
        monkeypatch.setitem(sys.modules, "pycsamt.emtools._core", dummy)

        dlg = FormatConverterDialog(small_sites)
        assert dlg._src_lbl.text() == "Loaded survey"
        dlg.close()


# ── FormatConverterDialog — _pick_dir ─────────────────────────────────────────


class TestPickDir:
    def test_pick_dir_updates_line_edit(self, qapp, monkeypatch, tmp_path):
        dlg = FormatConverterDialog(None)
        monkeypatch.setattr(
            converter_tool.QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: str(tmp_path)),
        )
        dlg._pick_dir()
        assert dlg._dir_edit.text() == str(tmp_path)
        dlg.close()

    def test_pick_dir_cancelled_leaves_line_edit_unchanged(
        self, qapp, monkeypatch
    ):
        dlg = FormatConverterDialog(None)
        original = dlg._dir_edit.text()
        monkeypatch.setattr(
            converter_tool.QFileDialog,
            "getExistingDirectory",
            staticmethod(lambda *a, **k: ""),
        )
        dlg._pick_dir()
        assert dlg._dir_edit.text() == original
        dlg.close()


# ── FormatConverterDialog — _on_run ───────────────────────────────────────────


class TestOnRunNoSites:
    def test_on_run_without_sites_logs_message_and_no_worker(self, qapp):
        dlg = FormatConverterDialog(None)
        dlg._on_run()
        assert "Load survey data first." in dlg._log.toPlainText()
        assert dlg._worker is None
        dlg.close()


class TestOnRunWithFakeWorker:
    def test_constructs_worker_with_current_ui_selections(
        self, qapp, monkeypatch, small_sites, tmp_path
    ):
        fake_cls = _fake_worker_cls(done_msg="Converted 3 stations")
        monkeypatch.setattr(converter_tool, "_ConvertWorker", fake_cls)

        dlg = FormatConverterDialog(small_sites)
        dlg._dir_edit.setText(str(tmp_path))
        dlg._fmt_combo.setCurrentText("JSON (station metadata)")
        dlg._on_run()

        assert len(fake_cls.captured) == 1
        worker = fake_cls.captured[0]
        assert worker.sites is small_sites
        assert worker.out_dir == tmp_path
        assert worker.fmt == "JSON (station metadata)"
        dlg.close()

    def test_run_button_disabled_during_run_reenabled_on_done(
        self, qapp, monkeypatch, small_sites, tmp_path
    ):
        states_during_progress = []

        def _capture_state(fake_self):
            states_during_progress.append(
                dlg._run_btn.isEnabled()
            )

        fake_cls = _fake_worker_cls(
            progress_events=[(1, 1, "S1")],
            done_msg="Converted 1 stations",
            on_start=_capture_state,
        )
        monkeypatch.setattr(converter_tool, "_ConvertWorker", fake_cls)

        dlg = FormatConverterDialog(small_sites)
        dlg._dir_edit.setText(str(tmp_path))
        assert dlg._run_btn.isEnabled()
        dlg._on_run()

        # Disabled at the moment the (fake) worker started running...
        assert states_during_progress == [False]
        # ...and re-enabled once `done` fired.
        assert dlg._run_btn.isEnabled()
        dlg.close()

    def test_log_cleared_and_progress_reset_on_run(
        self, qapp, monkeypatch, small_sites, tmp_path
    ):
        fake_cls = _fake_worker_cls(done_msg="Converted 3 stations")
        monkeypatch.setattr(converter_tool, "_ConvertWorker", fake_cls)

        dlg = FormatConverterDialog(small_sites)
        dlg._log.append("stale previous log line")
        dlg._dir_edit.setText(str(tmp_path))
        dlg._on_run()

        text = dlg._log.toPlainText()
        assert "stale previous log line" not in text
        assert f"Converting → {tmp_path}" in text
        dlg.close()

    def test_progress_events_update_bar_and_log(
        self, qapp, monkeypatch, small_sites, tmp_path
    ):
        fake_cls = _fake_worker_cls(
            progress_events=[(1, 2, "S1"), (2, 2, "S2")],
            done_msg="Converted 2 stations",
        )
        monkeypatch.setattr(converter_tool, "_ConvertWorker", fake_cls)

        dlg = FormatConverterDialog(small_sites)
        dlg._dir_edit.setText(str(tmp_path))
        dlg._on_run()

        assert dlg._progress.maximum() == 2
        assert dlg._progress.value() == 2
        text = dlg._log.toPlainText()
        assert "[1/2]" in text
        assert "S1" in text
        assert "[2/2]" in text
        assert "S2" in text
        dlg.close()

    def test_error_path_reenables_button_and_logs_error(
        self, qapp, monkeypatch, small_sites, tmp_path
    ):
        fake_cls = _fake_worker_cls(error_msg="boom: disk full")
        monkeypatch.setattr(converter_tool, "_ConvertWorker", fake_cls)

        dlg = FormatConverterDialog(small_sites)
        dlg._dir_edit.setText(str(tmp_path))
        dlg._on_run()

        assert dlg._run_btn.isEnabled()
        text = dlg._log.toPlainText()
        assert "Error: boom: disk full" in text
        dlg.close()


# ── FormatConverterDialog — direct slot unit tests ────────────────────────────


class TestOnProgressOnDoneOnError:
    def test_on_progress_sets_bar_and_appends_log(self, qapp):
        dlg = FormatConverterDialog(None)
        dlg._on_progress(4, 10, "kap109")
        assert dlg._progress.maximum() == 10
        assert dlg._progress.value() == 4
        assert "[4/10]" in dlg._log.toPlainText()
        assert "kap109" in dlg._log.toPlainText()
        dlg.close()

    def test_on_done_reenables_button_and_logs(self, qapp):
        dlg = FormatConverterDialog(None)
        dlg._run_btn.setEnabled(False)
        dlg._on_done("Converted 5 stations → /tmp/out")
        assert dlg._run_btn.isEnabled()
        assert "Done — Converted 5 stations" in dlg._log.toPlainText()
        dlg.close()

    def test_on_error_reenables_button_and_logs(self, qapp):
        dlg = FormatConverterDialog(None)
        dlg._run_btn.setEnabled(False)
        dlg._on_error("Cannot iterate survey: boom")
        assert dlg._run_btn.isEnabled()
        assert "Error: Cannot iterate survey: boom" in dlg._log.toPlainText()
        dlg.close()
