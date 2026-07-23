# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Extra tests for SectionPanel, on top of test_section_panel.py (which
covers construction, controls, and the basic set_result/clear cycle with
_MockResult). This file drives the load-from-disk flow, the topo overlay
branches, the compare-overlay contours, and the summary-fallback path.

pycsamt.topo.PYCSAMT_TOPO is a process-wide singleton; every test that
touches it snapshots/restores enabled+exaggeration+clip_below_surface
around itself so nothing leaks into other test files.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.panels.section_panel import SectionPanel
from pycsamt.topo.config import PYCSAMT_TOPO


class _MockMesh:
    x_nodes = [0, 500, 1000, 1500, 2000]
    z_nodes = [0, 100, 300, 600]
    station_x = [500, 1000, 1500]


class _MockResult:
    workdir = "/tmp/mock_inversion"
    n_iterations = 42
    final_rms = 0.987
    mesh = _MockMesh()
    summary = "MockResult: 42 iterations, RMS=0.987"
    rho_2d = np.log10(np.full((3, 4), 100.0))


class _MockResultNoSummary:
    """No .summary attribute -> exercises _show_summary's fallback branch."""

    workdir = "/tmp/ref_inversion"
    n_iterations = 10
    final_rms = 1.5
    mesh = _MockMesh()
    rho_2d = np.log10(np.full((3, 4), 50.0))


class _MockResultNoRMS:
    """No .final_rms -> exercises the RMS-annotation except branch."""

    workdir = "/tmp/no_rms"
    mesh = _MockMesh()
    summary = "no rms here"
    rho_2d = np.log10(np.full((3, 4), 100.0))


@pytest.fixture
def panel(qapp):
    p = SectionPanel()
    yield p
    p.close()


@pytest.fixture(autouse=True)
def restore_topo_singleton():
    saved = (
        PYCSAMT_TOPO.enabled,
        PYCSAMT_TOPO.exaggeration,
        PYCSAMT_TOPO.clip_below_surface,
        PYCSAMT_TOPO.interp_method,
    )
    yield
    (
        PYCSAMT_TOPO.enabled,
        PYCSAMT_TOPO.exaggeration,
        PYCSAMT_TOPO.clip_below_surface,
        PYCSAMT_TOPO.interp_method,
    ) = saved


# ── _on_load / _on_load_ref ───────────────────────────────────────────────


def test_on_load_cancelled_noop(panel, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: "")
    )
    panel._on_load()
    assert panel._result is None


def test_on_load_success_calls_load_from_workdir(panel, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: str(tmp_path)),
    )
    monkeypatch.setattr(
        "pycsamt.models.occam2d.InversionResult",
        lambda workdir: _MockResult(),
    )
    panel._on_load()
    assert panel._result is not None


def test_on_load_ref_cancelled_noop(panel, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: "")
    )
    panel._on_load_ref()
    assert panel._result_ref is None


def test_on_load_ref_success_enables_clear_compare(panel, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog,
        "getExistingDirectory",
        staticmethod(lambda *a, **k: str(tmp_path)),
    )
    monkeypatch.setattr(
        "pycsamt.models.occam2d.InversionResult",
        lambda workdir: _MockResultNoSummary(),
    )
    panel._on_load_ref()
    assert panel._result_ref is not None
    assert panel._btn_clear_compare.isEnabled()


def test_load_from_workdir_bad_result_shows_error(panel, tmp_path, monkeypatch):
    def _raise(workdir):
        raise RuntimeError("bad dir")

    monkeypatch.setattr(
        "pycsamt.models.occam2d.InversionResult", _raise
    )
    panel._load_from_workdir(str(tmp_path), is_reference=False)
    assert "Failed to load" in panel._summary.toPlainText()
    assert panel._result is None


def test_on_clear_ref_resets_state(panel):
    panel._result_ref = _MockResult()
    panel._btn_clear_compare.setEnabled(True)
    panel._on_clear_ref()
    assert panel._result_ref is None
    assert not panel._btn_clear_compare.isEnabled()


# ── result_loaded signal ──────────────────────────────────────────────────


def test_set_result_emits_result_loaded(panel):
    received = []
    panel.result_loaded.connect(received.append)
    panel.set_result(_MockResult())
    assert received == ["/tmp/mock_inversion"]


def test_set_result_emit_exception_is_swallowed(panel):
    class _BadWorkdir:
        def __str__(self):
            raise RuntimeError("boom")

    class _Result(_MockResult):
        workdir = _BadWorkdir()

    panel.set_result(_Result())  # must not raise


# ── _show_summary fallback ────────────────────────────────────────────────


def test_show_summary_fallback_when_no_summary_attr(panel):
    panel.set_result(_MockResultNoSummary())
    text = panel._summary.toPlainText()
    assert "Workdir" in text
    assert "/tmp/ref_inversion" in text


# ── _draw_section branches ────────────────────────────────────────────────


def test_draw_section_depth_clip(panel):
    panel._depth_max.setValue(200.0)
    panel.set_result(_MockResult())  # must not raise with clipping active


def test_draw_section_show_stations_unchecked(panel):
    panel._show_stations.setChecked(False)
    panel.set_result(_MockResult())  # must not raise


def test_draw_section_missing_rms_attr_swallowed(panel):
    panel.set_result(_MockResultNoRMS())  # must not raise; no title crash


def test_draw_section_missing_station_x_swallowed(panel):
    class _MeshNoStationX:
        x_nodes = _MockMesh.x_nodes
        z_nodes = _MockMesh.z_nodes

    class _Result(_MockResult):
        mesh = _MeshNoStationX()

    panel._show_stations.setChecked(True)
    panel.set_result(_Result())  # must not raise


# ── Compare overlay (contours) ────────────────────────────────────────────


def test_compare_overlay_draws_contours(panel):
    panel.set_result(_MockResult())
    panel._result_ref = _MockResultNoSummary()
    panel._redraw()  # must not raise


def test_compare_overlay_bad_ref_swallowed(panel):
    panel.set_result(_MockResult())
    panel._result_ref = object()  # missing rho_2d/mesh entirely
    panel._redraw()  # must not raise -> _draw_section_contours' except


# ── Topo toggle / exaggeration ────────────────────────────────────────────


def test_on_topo_toggled_enables_exag_spin(panel):
    panel._chk_topo.setChecked(True)
    assert panel._spin_exag.isEnabled()
    assert PYCSAMT_TOPO.enabled is True


def test_on_topo_toggled_off_disables_exag_spin(panel):
    panel._chk_topo.setChecked(True)
    panel._chk_topo.setChecked(False)
    assert not panel._spin_exag.isEnabled()
    assert PYCSAMT_TOPO.enabled is False


def test_on_topo_toggled_no_result_no_redraw_error(panel):
    panel._chk_topo.setChecked(True)  # _result is None -> just skips redraw


def test_on_exag_changed_updates_config(panel):
    panel._chk_topo.setChecked(True)
    panel._spin_exag.setValue(2.5)
    assert PYCSAMT_TOPO.exaggeration == pytest.approx(2.5)


def test_on_exag_changed_without_topo_checked_still_updates_config(panel):
    panel._spin_exag.setValue(3.0)
    assert PYCSAMT_TOPO.exaggeration == pytest.approx(3.0)


# ── Topo-draped rendering ──────────────────────────────────────────────────


def test_draw_section_topo_enabled_with_elevation(panel, monkeypatch):
    PYCSAMT_TOPO.enabled = True
    panel._chk_topo.setChecked(True)

    class _MeshWithElev(_MockMesh):
        station_elev = [100.0, 110.0, 105.0]
        station_names = ["S1", "S2", "S3"]

    class _Result(_MockResult):
        mesh = _MeshWithElev()

    monkeypatch.setattr(
        "pycsamt.topo.drape.interp_elev",
        lambda sx, elev, xc, method: np.zeros_like(xc),
    )

    def _fake_drape(xn_c, zn_km, rho_lin, elev, exaggeration, clip_above_surface):
        nz = len(zn_km)
        nx = len(xn_c)
        return xn_c, np.tile(zn_km[:, None], (1, nx)), rho_lin

    monkeypatch.setattr("pycsamt.topo.drape.drape_section", _fake_drape)
    monkeypatch.setattr(
        "pycsamt.topo.overlay.draw_topo_section", lambda *a, **k: None
    )

    panel.set_result(_Result())  # must not raise through the topo path


def test_draw_section_topo_enabled_missing_elevation_falls_back(panel):
    PYCSAMT_TOPO.enabled = True
    panel._chk_topo.setChecked(True)

    class _MeshNoStationX:
        x_nodes = _MockMesh.x_nodes
        z_nodes = _MockMesh.z_nodes

    class _Result(_MockResult):
        mesh = _MeshNoStationX()

    panel.set_result(_Result())  # station_x missing -> except -> flat fallback


# ── Export ────────────────────────────────────────────────────────────────


def test_on_export_opens_dialog_with_canvas_figure(panel, monkeypatch):
    captured = []

    class _FakeExport:
        def __init__(self, *a, **kw):
            captured.append(kw)

        def exec(self):
            return None

    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog", _FakeExport
    )
    panel._on_export()
    assert captured[0]["figure"] is panel._canvas.figure
