# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for MainWindow's action handlers: data loading, tool/dialog
launchers, settings profile I/O, recompute/correction/conversion bridges,
export, session persistence, and error paths.

``test_main_window.py`` covers construction and static structure; this
file drives the actual slots so the many locally-imported dialog/tool
classes get exercised too — each is monkeypatched at its origin module
(since MainWindow imports them lazily inside the handler body) with a
lightweight stand-in so no real (and possibly blocking) Qt dialog ever
opens under the offscreen platform.
"""

from __future__ import annotations

import enum
from types import SimpleNamespace
from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")


# ── Generic fake-dialog helpers ──────────────────────────────────────────────


class _DialogCode(enum.IntEnum):
    Rejected = 0
    Accepted = 1


def _make_dialog(exec_return=_DialogCode.Accepted, **extra_attrs):
    """Build a lightweight QDialog-shaped stand-in class.

    Every constructor call is recorded on the returned class's
    ``.captured`` list so tests can assert on how MainWindow invoked it,
    without pulling in the real (potentially heavy) dialog.
    """
    captured: list = []

    class _FakeSignal:
        def __init__(self):
            self.connected: list = []

        def connect(self, fn):
            self.connected.append(fn)

        def emit(self, *a, **k):
            for fn in list(self.connected):
                fn(*a, **k)

    class _Fake:
        DialogCode = _DialogCode

        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs
            for k, v in extra_attrs.items():
                setattr(self, k, _FakeSignal() if v == "signal" else v)
            captured.append(self)

        def exec(self):
            return exec_return

    _Fake.captured = captured
    return _Fake


# ── Window fixture ────────────────────────────────────────────────────────


@pytest.fixture
def window(qapp, monkeypatch, tmp_path):
    from pycsamt.app.desktop.models.session import SessionState

    monkeypatch.setattr(SessionState, "load", classmethod(lambda cls: cls()))
    # Redirect any real session.save() during tests away from the user's
    # actual ~/.pycsamt/session.json.
    monkeypatch.setattr(
        "pycsamt.app.desktop.models.session._SESSION_PATH",
        tmp_path / "session.json",
    )
    from pycsamt.app.desktop.main_window import MainWindow

    win = MainWindow()
    yield win
    win.close()


@pytest.fixture
def loaded_data(simulated_edi):
    """Real Sites + matching DataFrame via the non-Qt DataController."""
    from pycsamt.app.desktop.controllers.data_controller import (
        DataController,
    )

    dc = DataController()
    sites = dc.load([simulated_edi])
    return sites, dc.dataframe


# ── Data loading pipeline ────────────────────────────────────────────────────


def _fake_loader_cls(*, sites=None, df=None, error=None):
    class _Signal:
        def __init__(self):
            self._fn = None

        def connect(self, fn):
            self._fn = fn

        def emit(self, *a):
            if self._fn is not None:
                self._fn(*a)

    class _FakeLoader:
        def __init__(self, paths, parent=None):
            self._paths = list(paths)
            self.progress = _Signal()
            self.finished = _Signal()
            self.error = _Signal()
            self.data_controller = SimpleNamespace(dataframe=df)

        def start(self):
            self.progress.emit(50)
            if error is not None:
                self.error.emit(error)
            else:
                self.finished.emit(sites)

    return _FakeLoader


def test_open_files_rejected_dialog_does_nothing(window, monkeypatch):
    fake = _make_dialog(exec_return=_DialogCode.Rejected)
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.load_data_dlg.LoadDataDialog", fake
    )
    window._open_files()
    assert fake.captured  # constructed
    assert window._loader is None  # never proceeded to loading


def test_open_files_accepted_empty_paths_does_nothing(window, monkeypatch):
    fake = _make_dialog(exec_return=_DialogCode.Accepted, selected_paths=[])
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.load_data_dlg.LoadDataDialog", fake
    )
    window._open_files()
    assert window._loader is None


def test_open_files_accepted_triggers_load(
    window, monkeypatch, loaded_data, simulated_edi
):
    sites, df = loaded_data
    fake_dialog = _make_dialog(
        exec_return=_DialogCode.Accepted,
        selected_paths=[str(simulated_edi)],
    )
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.load_data_dlg.LoadDataDialog",
        fake_dialog,
    )
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.loader_worker.LoaderWorker",
        _fake_loader_cls(sites=sites, df=df),
    )
    window._open_files()
    assert window._session.last_data_dir == str(simulated_edi.parent)
    assert str(simulated_edi) in window._session.recent_files
    # _on_data_loaded ran synchronously via the fake loader's .start()
    assert "stations loaded" in window._status_file_lbl.text()


def test_start_loading_error_path(window, monkeypatch):
    monkeypatch.setattr(
        "pycsamt.app.desktop.workers.loader_worker.LoaderWorker",
        _fake_loader_cls(error="boom"),
    )
    window._start_loading(["/some/file.edi"])
    assert "ERROR: boom" in window._log_panel._text.toPlainText()
    assert window._status_ready_lbl.text() == "Error  ✕"


def test_on_data_loaded_with_real_dataframe(window, loaded_data):
    sites, df = loaded_data
    window._loader = SimpleNamespace(data_controller=SimpleNamespace(dataframe=df))
    # n_stations reads `self._controller.sites`, which `set_sites()` fills
    # in before firing this callback in the real pipeline — mirror that
    # rather than calling `_on_data_loaded` in isolation.
    window._controller.set_sites(sites)
    assert window._station_panel._table._model.rowCount() == len(df)
    assert f"{len(df)} stations loaded" == window._status_file_lbl.text()
    assert not window._progress_bar.isVisible()


def test_on_data_loaded_without_loader_only_updates_status(window):
    window._loader = None
    window._on_data_loaded([])
    assert window._status_file_lbl.text() == "0 stations loaded"


# ── Station selection ────────────────────────────────────────────────────────


def test_on_station_selected_updates_detail_card(window, loaded_data):
    sites, df = loaded_data
    window._controller.sites = sites
    station_id = df["ID"].iloc[0]
    window._on_station_selected(station_id)
    # No exception + detail card attempted an update (best-effort assertion:
    # the method ran to completion without raising).


def test_on_station_double_clicked_opens_profile(window, loaded_data):
    sites, df = loaded_data
    window._controller.sites = sites
    station_id = df["ID"].iloc[0]
    window._on_station_double_clicked(station_id)
    assert window._profile_win.isVisible()


def test_on_show_on_map_shows_map_window(window, loaded_data):
    sites, df = loaded_data
    station_id = df["ID"].iloc[0]
    window._on_show_on_map(station_id)
    assert window._map_win.isVisible()


# ── Panel window helpers ─────────────────────────────────────────────────────


def test_show_window_moves_and_shows(window):
    win = window._profile_win
    assert not win.isVisible()
    window._show_window(win)
    assert win.isVisible()


def test_panel_windows_lists_all_ten(window):
    assert len(window._panel_windows()) == 10


def test_open_agent_master_success(window, monkeypatch):
    result = SimpleNamespace(url="http://127.0.0.1:8765", started=True)
    monkeypatch.setattr(
        "pycsamt.app.desktop.main_window.launch_agent_master",
        lambda open_browser=True: result,
    )
    window._open_agent_master()
    assert "Starting Agent Master" in window._log_panel._text.toPlainText()


def test_open_agent_master_handles_exception(window, monkeypatch):
    def _raise(open_browser=True):
        raise RuntimeError("no port")

    monkeypatch.setattr("pycsamt.app.desktop.main_window.launch_agent_master", _raise)
    window._open_agent_master()
    assert "Could not launch Agent Master" in window._log_panel._text.toPlainText()


# ── Correction / conversion bridges ──────────────────────────────────────────


def test_on_corrections_committed(window, loaded_data):
    sites, _df = loaded_data
    window._on_corrections_committed(sites)
    assert "(corrected)" in window._status_file_lbl.text()


def test_on_conversion_committed_success(window, loaded_data):
    sites, _df = loaded_data
    window._on_conversion_committed(sites)
    assert "(converted)" in window._status_file_lbl.text()


def test_on_conversion_committed_bad_input_logs_error(window):
    class _Boom:
        def __iter__(self):
            raise RuntimeError("boom")

    window._on_conversion_committed(_Boom())
    text = window._log_panel._text.toPlainText()
    assert "Could not commit converted dataset" in text


# ── Forward -> Inversion bridge ──────────────────────────────────────────────


def test_on_forward_send_to_inversion(window):
    payload = {"dim": "1D", "resistivity": [100.0, 50.0]}
    window._on_forward_send_to_inversion(payload)
    assert window._inversion_win.isVisible()


def test_open_inversion_wizard_shows_window(window):
    window._open_inversion_wizard()
    assert window._inversion_win.isVisible()


def test_on_inversion_result_ready_no_result_noop(window):
    window._on_inversion_result_ready({"result": None})
    # no exception; nothing to assert beyond "ran cleanly"


def test_on_inversion_result_ready_with_result_dir(window):
    result = SimpleNamespace(result_dir=str(1) and "does/not/exist")
    window._on_inversion_result_ready({"result": result})
    text = window._log_panel._text.toPlainText()
    assert "Inversion model forwarded to Interpretation Studio." in text


# ── Export figure ────────────────────────────────────────────────────────────


def test_on_export_figure_no_figure_shows_status(window):
    window._on_export_figure()
    # nothing visible => early return; just must not raise


def test_on_export_figure_with_visible_canvas(window, monkeypatch):
    fake_figure = object()
    window._profile_win._canvas = SimpleNamespace(figure=fake_figure)
    window._profile_win.show()

    fake_export = _make_dialog(exec_return=_DialogCode.Accepted)
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog", fake_export
    )
    window._on_export_figure()
    assert fake_export.captured
    assert fake_export.captured[0].kwargs.get("figure") is fake_figure


# ── Preferences ───────────────────────────────────────────────────────────────


def test_open_preferences_accepted_applies_theme(window, monkeypatch):
    fake = _make_dialog(exec_return=_DialogCode.Accepted)
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.preferences_dlg.PreferencesDialog",
        fake,
    )
    window._open_preferences()
    assert "Preferences saved." in window._log_panel._text.toPlainText()


def test_open_preferences_rejected_does_nothing(window, monkeypatch):
    fake = _make_dialog(exec_return=_DialogCode.Rejected)
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.preferences_dlg.PreferencesDialog",
        fake,
    )
    window._open_preferences()
    assert "Preferences saved." not in window._log_panel._text.toPlainText()


# ── API configuration / settings ─────────────────────────────────────────────


def test_open_api_config(window, monkeypatch):
    fake = _make_dialog(exec_return=_DialogCode.Accepted, settings_changed="signal")
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.settings_dialog.APIConfigDialog", fake
    )
    # Avoid writing to the real ~/.pycsamt/settings.json default path.
    monkeypatch.setattr(window._settings_ctrl, "save", lambda *a, **k: None)
    window._open_api_config(tab="display")
    assert fake.captured
    assert fake.captured[0].kwargs.get("open_tab") == "display"


def test_on_settings_changed_logs_touched_keys(window):
    window._on_settings_changed(["station", "section", "view_controls"])
    assert (
        "Settings applied: station, section, view_controls"
        in window._log_panel._text.toPlainText()
    )


def test_reset_all_settings_confirmed(window):
    # no_modal_dialogs autouse fixture makes QMessageBox.question return Yes
    window._reset_all_settings()
    assert (
        "All API settings reset to package defaults."
        in window._log_panel._text.toPlainText()
    )


def test_reset_all_settings_declined(window, monkeypatch):
    from PySide6.QtWidgets import QMessageBox

    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(lambda *a, **k: QMessageBox.StandardButton.Cancel),
    )
    window._reset_all_settings()
    assert (
        "All API settings reset to package defaults."
        not in window._log_panel._text.toPlainText()
    )


def test_save_settings_profile(window, monkeypatch, tmp_path):
    target = tmp_path / "profile.json"
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *a, **k: (str(target), "")),
    )
    window._save_settings_profile()
    assert target.exists()
    assert "Settings profile saved" in window._log_panel._text.toPlainText()


def test_save_settings_profile_cancelled(window, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *a, **k: ("", "")),
    )
    window._save_settings_profile()
    assert "Settings profile saved" not in window._log_panel._text.toPlainText()


def test_load_settings_profile(window, monkeypatch, tmp_path):
    target = tmp_path / "profile.json"
    window._settings_ctrl.save(target)
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getOpenFileName",
        staticmethod(lambda *a, **k: (str(target), "")),
    )
    window._load_settings_profile()
    assert "Settings profile loaded" in window._log_panel._text.toPlainText()


def test_load_settings_profile_bad_file_logs_error(window, monkeypatch, tmp_path):
    target = tmp_path / "bad.json"
    target.write_text("not json")
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getOpenFileName",
        staticmethod(lambda *a, **k: (str(target), "")),
    )
    window._load_settings_profile()
    assert "Could not load settings profile" in window._log_panel._text.toPlainText()


def test_load_settings_profile_cancelled(window, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getOpenFileName",
        staticmethod(lambda *a, **k: ("", "")),
    )
    window._load_settings_profile()  # must not raise


# ── require_sites guard ──────────────────────────────────────────────────────


def test_require_sites_true_when_sites_loaded(window, loaded_data):
    sites, _df = loaded_data
    window._controller.sites = sites
    assert window._require_sites("X") is True


def test_require_sites_false_triggers_open_when_accepted(window, monkeypatch):
    window._controller.sites = None
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.no_data_dialog.NoDataDialog.require",
        classmethod(lambda cls, parent, tool_name="": True),
    )
    triggered = []
    monkeypatch.setattr(window._act_open, "trigger", lambda: triggered.append(1))
    assert window._require_sites("X") is False
    assert triggered == [1]


def test_require_sites_false_no_trigger_when_declined(window, monkeypatch):
    window._controller.sites = None
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.no_data_dialog.NoDataDialog.require",
        classmethod(lambda cls, parent, tool_name="": False),
    )
    triggered = []
    monkeypatch.setattr(window._act_open, "trigger", lambda: triggered.append(1))
    assert window._require_sites("X") is False
    assert triggered == []


# ── Tool menu launchers ───────────────────────────────────────────────────────

_GATED_TOOLS = [
    (
        "_open_strike_analyzer",
        "pycsamt.app.desktop.tools.strike_tool.StrikeAnalyzerDialog",
    ),
    (
        "_open_format_converter",
        "pycsamt.app.desktop.tools.converter_tool.FormatConverterDialog",
    ),
    (
        "_open_station_response",
        "pycsamt.app.desktop.tools.station_response_tool.StationResponseDialog",
    ),
    (
        "_open_strike_profile",
        "pycsamt.app.desktop.tools.strike_profile_tool.StrikeProfileDialog",
    ),
    (
        "_open_phase_tensor_map",
        "pycsamt.app.desktop.tools.phase_tensor_map_tool.PhaseTensorMapDialog",
    ),
    (
        "_open_phase_tensor_strip_grid",
        "pycsamt.app.desktop.tools.phase_tensor_strip_grid_tool."
        "PhaseTensorStripGridDialog",
    ),
    (
        "_open_dimensionality",
        "pycsamt.app.desktop.tools.dimensionality_tool.DimensionalityDialog",
    ),
    (
        "_open_elevation_enrichment",
        "pycsamt.app.desktop.tools.elevation_tool.ElevationEnrichDialog",
    ),
]

_UNGATED_TOOLS = [
    (
        "_open_batch_export",
        "pycsamt.app.desktop.tools.batch_export_tool.BatchExportDialog",
    ),
    (
        "_open_coord_transformer",
        "pycsamt.app.desktop.tools.coord_tool.CoordTransformDialog",
    ),
    (
        "_open_layered_model",
        "pycsamt.app.desktop.tools.layered_model_tool.LayeredModelDialog",
    ),
]


@pytest.mark.parametrize("method_name,dotted_path", _GATED_TOOLS)
def test_gated_tool_blocked_without_sites(
    window, monkeypatch, method_name, dotted_path
):
    window._controller.sites = None
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.no_data_dialog.NoDataDialog.require",
        classmethod(lambda cls, parent, tool_name="": False),
    )
    fake = _make_dialog()
    monkeypatch.setattr(dotted_path, fake)
    getattr(window, method_name)()
    assert fake.captured == []


@pytest.mark.parametrize("method_name,dotted_path", _GATED_TOOLS)
def test_gated_tool_opens_with_sites(
    window, monkeypatch, loaded_data, method_name, dotted_path
):
    sites, _df = loaded_data
    window._controller.sites = sites
    fake = _make_dialog()
    monkeypatch.setattr(dotted_path, fake)
    getattr(window, method_name)()
    assert len(fake.captured) == 1


@pytest.mark.parametrize("method_name,dotted_path", _UNGATED_TOOLS)
def test_ungated_tool_always_opens(window, monkeypatch, method_name, dotted_path):
    fake = _make_dialog()
    monkeypatch.setattr(dotted_path, fake)
    getattr(window, method_name)()
    assert len(fake.captured) == 1


def test_edi_validator_applies_modified_sites(window, monkeypatch, loaded_data):
    sites, _df = loaded_data
    window._controller.sites = sites
    fake = _make_dialog(
        exec_return=1,
        modified_sites=sites,
        open_recompute_requested="signal",
    )
    monkeypatch.setattr(
        "pycsamt.app.desktop.tools.validator_tool.EDIValidatorDialog", fake
    )
    window._open_edi_validator()
    assert "EDI Validator" in window._log_panel._text.toPlainText()


def test_edi_validator_no_modification_noop(window, monkeypatch, loaded_data):
    sites, _df = loaded_data
    window._controller.sites = sites
    fake = _make_dialog(
        exec_return=0, modified_sites=None, open_recompute_requested="signal"
    )
    monkeypatch.setattr(
        "pycsamt.app.desktop.tools.validator_tool.EDIValidatorDialog", fake
    )
    window._open_edi_validator()
    assert fake.captured


def test_frequency_editor_applies_edited_sites(window, monkeypatch, loaded_data):
    sites, _df = loaded_data
    window._controller.sites = sites
    fake = _make_dialog(exec_return=1, edited_sites=sites)
    monkeypatch.setattr(
        "pycsamt.app.desktop.tools.frequency_editor_tool.FrequencyEditorDialog",
        fake,
    )
    window._open_frequency_editor()
    assert "Frequency Editor" in window._log_panel._text.toPlainText()


def test_collect_figures_empty_by_default(window):
    assert window._collect_figures() == []


def test_collect_figures_finds_visible_canvas(window):
    window._profile_win._canvas = SimpleNamespace(figure=object())
    window._profile_win.show()
    figs = window._collect_figures()
    assert any(label == "Profile" for label, _fig in figs)


# ── Recompute ─────────────────────────────────────────────────────────────────


def test_open_recompute_no_prior_recompute(window, monkeypatch):
    fake = _make_dialog(recompute_committed="signal")
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.recompute_dlg.RecomputeDialog", fake
    )
    window._open_recompute()
    assert fake.captured


def test_open_recompute_confirmed_when_already_recomputed(window, monkeypatch):
    window._recomputed_ids = {"S1"}
    fake = _make_dialog(recompute_committed="signal")
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.recompute_dlg.RecomputeDialog", fake
    )
    # no_modal_dialogs fixture returns Yes by default -> proceeds
    window._open_recompute()
    assert fake.captured


def test_open_recompute_declined_when_already_recomputed(window, monkeypatch):
    window._recomputed_ids = {"S1"}
    from PySide6.QtWidgets import QMessageBox

    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(lambda *a, **k: QMessageBox.StandardButton.No),
    )
    fake = _make_dialog(recompute_committed="signal")
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.recompute_dlg.RecomputeDialog", fake
    )
    window._open_recompute()
    assert fake.captured == []


def test_on_recompute_committed(window, loaded_data, tmp_path):
    sites, _df = loaded_data
    result = SimpleNamespace(
        records=[
            SimpleNamespace(station="S1", status="ok"),
            SimpleNamespace(station="S2", status="failed"),
        ],
        output_root=tmp_path,
        sites=sites,
    )
    window._on_recompute_committed(result)
    # NOTE (found while writing this test, not fixed here — flagged
    # separately): `_apply_modified_sites` re-triggers `_on_data_loaded`
    # via `_controller.set_sites`, which unconditionally clears
    # `_recomputed_ids` *after* this method has already populated it and
    # *before* `mark_recomputed` reads it — so the badge set is empty by
    # the time it reaches the table model. Asserting actual behavior here.
    assert window._recomputed_ids == set()
    assert window._last_recompute_output == tmp_path
    assert "1 station(s) marked with" in window._log_panel._text.toPlainText()
    assert "marked with" in window._log_panel._text.toPlainText()


# ── Help menu ─────────────────────────────────────────────────────────────────


def test_open_documentation_opens_url(window, monkeypatch):
    from PySide6.QtGui import QDesktopServices

    calls = []
    monkeypatch.setattr(
        QDesktopServices,
        "openUrl",
        staticmethod(lambda url: calls.append(url.toString())),
    )
    window._open_documentation()
    assert calls == ["https://pycsamt.org/"]


def test_open_github_opens_url(window, monkeypatch):
    from PySide6.QtGui import QDesktopServices

    calls = []
    monkeypatch.setattr(
        QDesktopServices,
        "openUrl",
        staticmethod(lambda url: calls.append(url.toString())),
    )
    window._open_github()
    assert calls == ["https://github.com/earthai-tech/pycsamt"]


def test_open_about(window, monkeypatch):
    fake = _make_dialog()
    monkeypatch.setattr("pycsamt.app.desktop.dialogs.about_dialog.AboutDialog", fake)
    window._open_about()
    assert fake.captured


# ── Recent files ──────────────────────────────────────────────────────────────


def test_rebuild_recent_menu_empty(window):
    window._session.recent_files = []
    window._rebuild_recent_menu()
    actions = window._act_recent.actions()
    assert len(actions) == 1
    assert actions[0].text() == "(none)"
    assert not actions[0].isEnabled()


def test_rebuild_recent_menu_triggers_load(window, monkeypatch):
    window._session.recent_files = ["/data/A.edi", "/data/B.edi"]
    window._rebuild_recent_menu()
    actions = window._act_recent.actions()
    assert [a.text() for a in actions] == ["/data/A.edi", "/data/B.edi"]

    called = []
    monkeypatch.setattr(window, "_start_loading", lambda paths: called.append(paths))
    actions[0].trigger()
    assert called == [["/data/A.edi"]]


# ── Filter / search ───────────────────────────────────────────────────────────


def test_on_filter_changed_delegates(window, monkeypatch):
    # StationPanel has no `.filter()` method today, so the `try` line
    # itself would raise AttributeError without this stand-in — add one
    # dynamically (raising=False: the attribute doesn't exist yet) to
    # exercise the success path.
    calls = []
    monkeypatch.setattr(
        window._station_panel,
        "filter",
        lambda text: calls.append(text),
        raising=False,
    )
    window._on_filter_changed("WILLY")
    assert calls == ["WILLY"]


def test_on_filter_changed_swallows_exception(window):
    # StationPanel currently has no `.filter()` method, so this always
    # takes the except branch in production too — just confirms it never
    # propagates.
    window._on_filter_changed("x")  # must not raise


# ── Session ───────────────────────────────────────────────────────────────────


def test_on_save_session(window):
    window._on_save_session()
    assert "Session saved." in window._log_panel._text.toPlainText()


def test_restore_layout_handles_bad_geometry(window):
    window._session.dock_geometry = 12345  # not a str -> .encode() raises
    window._session.dock_state = 12345
    window._restore_layout()  # must not raise


# ── Error handling ────────────────────────────────────────────────────────────


def test_on_load_error_updates_ui(window):
    window._on_load_error("disk full")
    assert not window._progress_bar.isVisible()
    assert window._status_ready_lbl.text() == "Error  ✕"
    assert "ERROR: disk full" in window._log_panel._text.toPlainText()


# ── closeEvent ────────────────────────────────────────────────────────────────


def test_close_event_declined_ignores_and_keeps_open(window, monkeypatch):
    from PySide6.QtWidgets import QMessageBox

    monkeypatch.setattr(
        QMessageBox,
        "question",
        staticmethod(lambda *a, **k: QMessageBox.StandardButton.Cancel),
    )
    saved = []
    monkeypatch.setattr(window._session, "save", lambda: saved.append(1))
    event = mock.Mock()
    window.closeEvent(event)
    event.ignore.assert_called_once()
    assert saved == []


# ── Remaining coverage gaps ──────────────────────────────────────────────────


def test_edi_validator_blocked_without_sites(window, monkeypatch):
    window._controller.sites = None
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.no_data_dialog.NoDataDialog.require",
        classmethod(lambda cls, parent, tool_name="": False),
    )
    fake = _make_dialog()
    monkeypatch.setattr(
        "pycsamt.app.desktop.tools.validator_tool.EDIValidatorDialog", fake
    )
    window._open_edi_validator()
    assert fake.captured == []


def test_frequency_editor_blocked_without_sites(window, monkeypatch):
    window._controller.sites = None
    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.no_data_dialog.NoDataDialog.require",
        classmethod(lambda cls, parent, tool_name="": False),
    )
    fake = _make_dialog()
    monkeypatch.setattr(
        "pycsamt.app.desktop.tools.frequency_editor_tool.FrequencyEditorDialog",
        fake,
    )
    window._open_frequency_editor()
    assert fake.captured == []


def test_on_inversion_result_ready_with_rho_2d_model(window):
    """Exercises the ``hasattr(result, "rho_2d")`` branch (as opposed to
    the ``result_dir`` branch already covered elsewhere). Attributes
    below are all populated deliberately: a model missing e.g.
    ``depth_max`` trips the unrelated, already-documented
    ``_update_status_card`` crash bug (see
    test_load_from_inversion_model_missing_depth_max_crashes in
    test_interp_window.py)."""
    model = SimpleNamespace(
        rho_2d=[[1.0]],
        n_x=4,
        n_z=3,
        depth_max=100.0,
        profile_length=500.0,
        method="Occam2D",
        rms=1.0,
    )
    window._on_inversion_result_ready({"result": model})
    assert window._interp_win._ctrl.state.model is model
    text = window._log_panel._text.toPlainText()
    assert "Inversion model forwarded to Interpretation Studio." in text


def test_on_data_loaded_panel_window_failure_does_not_block_others(
    window, loaded_data, monkeypatch
):
    """One panel-window setter raising must not stop the others from
    being fed, nor block the final status-bar update."""
    sites, df = loaded_data
    window._loader = SimpleNamespace(data_controller=SimpleNamespace(dataframe=df))

    def _boom(_sites):
        raise RuntimeError("panel boom")

    monkeypatch.setattr(window._profile_win, "set_sites", _boom)
    calls = []
    monkeypatch.setattr(window._map_win, "set_sites", lambda s: calls.append(s))
    # set_sites() on the controller is what actually fires _on_data_loaded
    # (wired via ctrl.on_data_loaded in __init__); calling the handler
    # directly would leave ctrl.n_stations at 0 for the final assertion.
    window._controller.set_sites(sites)
    assert calls == [sites]
    assert f"{len(df)} stations loaded" in window._status_file_lbl.text()


def test_apply_theme_panel_window_failure_swallowed(window, monkeypatch):
    def _boom(_dark):
        raise RuntimeError("theme boom")

    monkeypatch.setattr(window._profile_win, "set_dark_mode", _boom)
    window._apply_theme("dark")  # must not raise
    assert window._session.theme == "dark"


def test_on_station_selected_propagates_to_visible_windows(window, loaded_data):
    sites, df = loaded_data
    window._controller.sites = sites
    station_id = df["ID"].iloc[0]
    window._profile_win.show()
    window._map_win.show()
    window._on_station_selected(station_id)  # must not raise


def test_on_settings_changed_refreshes_visible_profile_and_qc_windows(
    window, monkeypatch
):
    window._profile_win.show()
    window._qc_win.show()
    profile_calls = []
    qc_calls = []
    monkeypatch.setattr(
        window._profile_win, "_on_refresh", lambda: profile_calls.append(1)
    )
    monkeypatch.setattr(window._qc_win, "_on_run", lambda: qc_calls.append(1))
    window._on_settings_changed(["station", "section"])
    assert profile_calls == [1]
    assert qc_calls == [1]


def test_on_settings_changed_swallows_refresh_exceptions(window):
    window._profile_win.show()
    window._qc_win.show()
    window._profile_win._on_refresh = mock.Mock(side_effect=RuntimeError("boom"))
    window._qc_win._on_run = mock.Mock(side_effect=RuntimeError("boom"))
    window._on_settings_changed(["station", "section"])  # must not raise


def test_save_settings_profile_error_logged(window, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *a, **k: ("/fake/out.json", "")),
    )

    def _boom(path):
        raise RuntimeError("disk full")

    monkeypatch.setattr(window._settings_ctrl, "save", _boom)
    window._save_settings_profile()
    assert "Save profile error: disk full" in window._log_panel._text.toPlainText()
