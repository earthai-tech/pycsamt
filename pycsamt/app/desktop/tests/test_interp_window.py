# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the Interpretation Studio window
(pycsamt.app.desktop.windows.interp).

Most tests use a small synthetic ResistivityModel so they stay fast; the
real-data tests load the bundled Occam2D run folder and a ModEM 3-D PCSF
through the window's background loader, exactly as a user would.
QFileDialog / QInputDialog are never reached: every evidence action takes
its path (and PCBH positions) as arguments.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers.interp_studio import (
    NEEDS,
    TABS,
    tab,
    view_status,
)
from pycsamt.app.desktop.windows.interp import InterpretationWindow

ROOT = Path(__file__).resolve().parents[4]
OCCAM2D = ROOT / "data" / "occam2D"
MODEM3D_PCSF = (ROOT / "examples" / "pcsf_conversion_demo" / "output"
                / "modem3d_no_topo.pcsf")
PCBH = (ROOT / "pycsamt" / "format" / "tests" / "data" / "borehole"
        / "minimal_vertical.pcbh.json")


def _model(names=("S1", "S2", "S3")):
    from pycsamt.interp import ResistivityModel

    rho = np.log10(np.array([[800.0, 800.0, 800.0],
                             [600.0, 800.0, 700.0],
                             [15.0, 800.0, 20.0],
                             [10.0, 800.0, 12.0],
                             [8.0, 800.0, 9.0]]))
    return ResistivityModel.from_array(
        rho, x_centers=np.array([0.0, 500.0, 1000.0]),
        z_centers=np.array([5.0, 15.0, 30.0, 60.0, 100.0]),
        station_x=np.array([0.0, 500.0, 1000.0]),
        station_names=list(names), method="TDEM", rms=1.5)


@pytest.fixture
def win(qapp):
    w = InterpretationWindow(parent=None)
    yield w
    w.wait_for_task()
    w.close()
    import matplotlib.pyplot as plt

    plt.close("all")


@pytest.fixture
def loaded(win):
    win.set_model(_model(), "synthetic")
    return win


def _shown(w) -> bool:
    return w._view._stack.currentWidget() is w._view.canvas


def _card_text(w) -> str:
    from PySide6.QtWidgets import QLabel

    return " ".join(lbl.text() for lbl in
                    w._view._unavailable.findChildren(QLabel))


# ── catalogue ───────────────────────────────────────────────────────────────


class TestCatalogue:
    def test_every_view_is_a_controller_method(self):
        from pycsamt.app.desktop.controllers.interp_controller import (
            InterpController,
        )

        for t in TABS:
            for v in t.views:
                assert callable(getattr(InterpController, v.method, None)), \
                    v.method
                assert set(v.needs) <= set(NEEDS)

    def test_steps_are_controller_methods(self):
        from pycsamt.app.desktop.controllers.interp_controller import (
            InterpController,
        )

        for t in TABS:
            if t.step:
                assert callable(getattr(InterpController, t.step))

    def test_view_status_reports_what_is_missing(self):
        from pycsamt.app.desktop.controllers.interp_controller import (
            InterpController,
        )

        c = InterpController()
        v = tab("hydrology").views[0]
        ok, why = view_status(c, v)
        assert not ok and "hydrology" in why.lower()


# ── construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_title_and_tabs(self, win):
        assert "Interpretation Studio" in win.windowTitle()
        assert [b.text() for b in win._tab_buttons] == [t.label
                                                         for t in TABS]
        assert win._tab_key == "geology"

    def test_empty_state_explains_itself(self, win):
        assert not _shown(win)
        assert "model" in _card_text(win).lower()
        assert "No model" in win._evidence.model_card.text()

    def test_figures_are_publication_white(self, win):
        win.set_dark_mode(True)
        assert win._ctrl.dark is False

    def test_settings_forms_only_for_tabs_with_fields(self, win):
        assert set(win._forms) == {t.key for t in TABS if t.fields}

    def test_views_marked_missing_without_model(self, win):
        from pycsamt.app.desktop.windows.interp.window import _MISSING

        item = win._view_list.item(0)
        assert item.foreground().color() == _MISSING
        assert "Needs:" in item.toolTip()

    def test_step_button_label_keeps_ampersand(self, win):
        # a single "&" would be eaten as a Qt mnemonic
        assert "&&" in win._btn_step.text()
        win._select_tab(1)  # Structure has no step
        assert not win._btn_step.isVisibleTo(win)


# ── model + views ───────────────────────────────────────────────────────────


class TestModelAndViews:
    def test_model_view_draws(self, loaded):
        assert _shown(loaded)
        assert "3 × 5" in loaded._evidence.model_card.text()
        assert loaded._station_combo.count() == 3

    def test_needs_run_view_shows_reason_card(self, loaded):
        loaded.select_view("plot_strat_log")
        assert loaded.current_view() == "plot_strat_log"
        assert not _shown(loaded)
        assert "classify" in _card_text(loaded).lower()

    def test_geology_step_unlocks_logs(self, loaded):
        loaded._run_step()
        assert loaded.wait_for_task()
        assert len(loaded._ctrl.state.strat_logs) == 3
        assert "Classified 3" in loaded._step_status.text()
        loaded.select_view("plot_strat_log")
        assert _shown(loaded)

    def test_hydrology_step_uses_form_settings(self, loaded):
        loaded._forms["hydrology"].set_values(
            {"petro_model": "waxman_smits", "rho_w": 25.0, "sigma_s": 0.2})
        loaded._select_tab(2)
        loaded._run_step()
        assert loaded.wait_for_task()
        cfg = loaded._ctrl.state.petro_cfg
        assert type(cfg.petro).__name__ == "WaxmanSmitsModel"
        assert cfg.rho_w == pytest.approx(25.0)
        assert loaded._ctrl.state.hydro_result is not None
        assert _shown(loaded)  # K section

    def test_step_without_model_is_refused(self, win):
        win._run_step()
        assert win._task is None
        assert "model" in win._step_status.text().lower()

    def test_per_station_view_follows_station(self, loaded):
        loaded._run_step()
        loaded.wait_for_task()
        loaded.select_view("plot_strat_log")
        loaded._station_combo.setCurrentIndex(2)
        loaded._on_station(2)
        title = loaded._figure.texts[0].get_text() if \
            loaded._figure.texts else loaded._figure._suptitle.get_text()
        assert "S3" in title

    def test_tab_remembers_its_view(self, loaded):
        loaded.select_view("plot_rock_db")
        loaded._select_tab(5)
        loaded._select_tab(0)
        assert loaded.current_view() == "plot_rock_db"

    def test_superseded_figures_are_closed(self, loaded):
        import matplotlib.pyplot as plt

        for _ in range(4):
            loaded._rerender()
        assert len(plt.get_fignums()) <= 2


# ── evidence ────────────────────────────────────────────────────────────────


class TestEvidence:
    def test_pcbh_borehole_at_matching_station(self, qapp):
        w = InterpretationWindow()
        from pycsamt.format.borehole import read_pcbh

        bid = read_pcbh(PCBH).boreholes[0].id
        w.set_model(_model(names=("S1", bid, "S3")))
        w._evidence.add_borehole_pcbh(str(PCBH))  # no dialog: name matches
        bh = w._ctrl.state.boreholes
        assert [b.name for b in bh] == [bid] or len(bh) == 1
        assert bh[0].x == pytest.approx(500.0)
        assert w._evidence.bh_list.count() == 1
        w.close()

    def test_pcbh_explicit_position(self, loaded):
        from pycsamt.format.borehole import read_pcbh

        bid = read_pcbh(PCBH).boreholes[0].id
        loaded._evidence.add_borehole_pcbh(str(PCBH),
                                           positions={bid: 250.0})
        assert loaded._ctrl.state.boreholes[0].x == pytest.approx(250.0)
        ok, _ = view_status(loaded._ctrl, tab("geology").views[6])
        assert ok  # borehole fence is ready

    def test_remove_borehole(self, loaded):
        from pycsamt.format.borehole import read_pcbh

        bid = read_pcbh(PCBH).boreholes[0].id
        loaded._evidence.add_borehole_pcbh(str(PCBH), positions={bid: 0.0})
        loaded._evidence.bh_list.setCurrentRow(0)
        loaded._evidence.remove_borehole()
        assert not loaded._ctrl.state.boreholes

    def test_bad_file_is_reported_not_raised(self, loaded, tmp_path):
        bad = tmp_path / "x.csv"
        bad.write_text("nonsense")
        msgs = []
        loaded._evidence.changed.connect(msgs.append)
        loaded._evidence.add_borehole_csv(str(bad))
        assert msgs and msgs[-1].startswith("Failed")

    def test_repeat_survey_needs_baseline(self, win):
        msgs = []
        win._evidence.changed.connect(msgs.append)
        win._evidence.add_survey(str(OCCAM2D), label="2025")
        assert "baseline" in msgs[-1].lower()


# ── library & export ────────────────────────────────────────────────────────


class TestLibrary:
    def test_pin_and_export(self, loaded, tmp_path):
        loaded._pin_current()
        loaded.select_view("plot_rock_db")
        loaded._pin_current()
        assert len(loaded.pinned_labels()) == 2
        assert loaded._gallery.count() == 2
        paths = loaded.export_pinned(tmp_path, fmt="png", dpi=60)
        assert len(paths) == 2 and all(p.stat().st_size > 0 for p in paths)

    def test_pinned_figure_survives_view_change(self, loaded):
        fig = loaded._figure
        loaded._pin_current()
        loaded.select_view("plot_rock_db")
        loaded._show_pinned(loaded._gallery.item(0))
        assert loaded._figure is fig

    def test_unpin(self, loaded):
        loaded._pin_current()
        loaded._gallery.setCurrentRow(0)
        loaded._unpin()
        assert loaded._gallery.count() == 0 and not loaded.pinned_labels()


# ── session + main-window API ───────────────────────────────────────────────


class TestSession:
    def test_settings_and_tab_round_trip(self, win, qapp):
        win._forms["uncertainty"].set_values({"n_samples": 77})
        win._select_tab(4)
        store = {}
        win.save_geometry_to(store)
        w2 = InterpretationWindow()
        w2.restore_geometry_from(store)
        assert w2._tab_key == "uncertainty"
        assert w2._forms["uncertainty"].values()["n_samples"] == 77
        w2.close()

    def test_close_hides_and_signals(self, win):
        got = []
        win.panel_closed.connect(lambda: got.append(1))
        win.show()
        win.close()
        assert got and not win.isVisible()

    def test_receive_inversion_in_memory_model(self, win):
        win.receive_inversion({"engine": "x", "result": _model()})
        assert win._ctrl.state.model is not None
        assert _shown(win)

    def test_legacy_import_path(self):
        from pycsamt.app.desktop.windows import interp_window

        assert interp_window.InterpretationWindow is InterpretationWindow


# ── real data ───────────────────────────────────────────────────────────────


@pytest.mark.skipif(not OCCAM2D.is_dir(), reason="Occam2D sample missing")
class TestRealData:
    def test_occam2d_run_folder_through_inversion_payload(self, win):
        win.receive_inversion({"engine": "occam2d", "result": None,
                               "path": str(OCCAM2D)})
        assert win.wait_for_task()
        m = win._ctrl.state.model
        assert m is not None and m.n_x > 100
        assert "Occam2D" in win._evidence.model_card.text()
        assert _shown(win)
        win._run_step()
        assert win.wait_for_task()
        assert len(win._ctrl.state.strat_logs) == len(m.station_names)

    def test_repeat_survey_on_same_grid(self, win):
        win.load_model(str(OCCAM2D))
        win.wait_for_task()
        win._evidence.add_survey(str(OCCAM2D), label="2025")
        assert win._ctrl.state.timelapse_labels == ["baseline", "2025"]
        win.select_view("plot_timelapse_change")
        assert _shown(win)

    @pytest.mark.skipif(not MODEM3D_PCSF.is_file(), reason="no PCSF sample")
    def test_modem3d_pcsf_lines(self, win):
        win.load_model(str(MODEM3D_PCSF))
        assert win.wait_for_task()
        info = win._ctrl.state.model_info
        assert info is not None and info.kind == "pcsf"
        many = len(info.lines) > 1
        assert win._evidence.line_combo.isVisibleTo(win._evidence) == many
        if many:
            win.load_model(str(MODEM3D_PCSF), info.lines[1])
            assert win.wait_for_task()
            assert win._ctrl.state.model_info.line == info.lines[1]

    def test_bad_source_is_reported(self, win, tmp_path):
        win.load_model(str(tmp_path))
        assert win.wait_for_task()
        assert win._ctrl.state.model is None
        assert "not loaded" in win._step_status.text().lower()
