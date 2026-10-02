# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for Pcsf3DWindow (pycsamt.app.desktop.windows.pcsf3d_window).

The desktop app's 3-D/volume panel: loads a ``.pcsf`` file via
``pycsamt.map.MapView.from_pcsf`` and renders it through
``PlotlyView`` (QtWebEngine). Loads a real Occam2D-derived .pcsf file
(built from the real bundled data/occam2D dataset) so the rendered scene
is genuine inverted resistivity, not synthetic placeholder data.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")
pytest.importorskip(
    "PySide6.QtWebEngineWidgets", reason="QtWebEngine required"
)

from PySide6.QtWidgets import QWidget

from pycsamt.app.desktop.windows.pcsf3d_window import Pcsf3DWindow

_ROOT = Path(__file__).parents[4]
_OCCAM_DIR = _ROOT / "data" / "occam2D"
_SKIP_OCCAM = pytest.mark.skipif(
    not _OCCAM_DIR.exists(), reason=f"bundled occam2D data not found: {_OCCAM_DIR}"
)


def _write_grid2d_pcsf(tmp_path) -> Path:
    from pycsamt.format import write_pcsf
    from pycsamt.format.adapters.occam2d import occam2d_to_pcsf
    from pycsamt.models.occam2d.results import InversionResult

    result = InversionResult(workdir=_OCCAM_DIR)
    model = occam2d_to_pcsf(result)
    return write_pcsf(model, tmp_path / "grid2d.pcsf")


def _geo_referenced_view():
    """A MapView with real lat/lon stations on two lines -- mirrors
    pycsamt/map/tests/test_borehole_align.py's own survey fixture. No
    resistivity/profile data, so map3d() alone renders no data traces;
    this is only meant to exercise the borehole-alignment plumbing, not
    a real fence render (the occam2D-derived fixture above has no real
    lat/lon of its own -- an Occam2D result is profile-only unless
    geo-referenced against a real EDI survey via known_stations)."""
    from pycsamt.map import MapView
    from pycsamt.map._core import MapData, StationRecord

    stations = []
    idx = 0
    for line_idx, lat in enumerate((22.000, 22.002)):
        for j in range(6):
            stations.append(
                StationRecord(
                    id=f"L{line_idx}-{j}",
                    latitude=lat,
                    longitude=103.000 + 0.001 * j,
                    elevation=500.0,
                    line=f"L{line_idx}",
                    index=idx,
                )
            )
            idx += 1
    return MapView(MapData(sites=None, stations=tuple(stations)))


def _pcbh_document():
    from pycsamt.format.borehole import (
        Collar,
        CoordinateReferenceSystem,
        LogInterval,
        PCBHBorehole,
        PCBHDocument,
        VocabularyEntry,
    )

    collar = Collar(
        103.001, 22.0005, 500.0, longitude=103.001, latitude=22.0005
    )
    return PCBHDocument(
        document_id="test:pcsf3d",
        created_at="2026-09-22T00:00:00Z",
        created_by="pytest",
        crs=CoordinateReferenceSystem("EPSG:4326"),
        lithologies=[
            VocabularyEntry("A", "Overburden", color="#cccccc"),
            VocabularyEntry("B", "Granite", color="#884422"),
        ],
        boreholes=[
            PCBHBorehole(
                id="BH1",
                name="Borehole 1",
                kind="mining_exploration",
                status="completed",
                collar=collar,
                total_depth_md=120.0,
                interval_logs={
                    "lithology": [
                        LogInterval(0.0, 40.0, code="A"),
                        LogInterval(40.0, 120.0, code="B"),
                    ]
                },
            )
        ],
    )


def _write_pcbh(tmp_path) -> Path:
    from pycsamt.format.borehole import write_pcbh

    return write_pcbh(_pcbh_document(), tmp_path / "boreholes.pcbh.json")


def _pcgs_document():
    from pycsamt.format.structure import StructModel
    from pycsamt.geology.structural import (
        FaultTrace,
        LinearMeasurement,
        StructuralMeasurement,
        StructuralModel,
    )

    model = StructuralModel(
        faults=[
            FaultTrace(
                x=50.0, dip_deg=60.0, downthrown_side="right",
                sense="normal", line="L0",
            )
        ],
        planar=[
            StructuralMeasurement(
                x=30.0, kind="bedding", strike_deg=45.0, dip_deg=30.0,
                dip_direction_deg=135.0, line="L0",
            )
        ],
        linear=[
            LinearMeasurement(
                x=70.0, kind="fold_axis", trend_deg=90.0, plunge_deg=20.0,
                line="L1",
            )
        ],
    )
    return StructModel(
        document_id="test:pcsf3d-struct",
        created_at="2026-09-22T00:00:00Z",
        created_by="pytest",
        model=model,
    )


def _write_pcgs(tmp_path) -> Path:
    from pycsamt.format.structure import write_structure

    return write_structure(_pcgs_document(), tmp_path / "structure.pcgs.json")


def _write_legend_pcgl(tmp_path) -> Path:
    from pycsamt.format import GeologyLegend, write_legend
    from pycsamt.geology import RockDatabase, RockEntry

    db = RockDatabase(
        [
            RockEntry("Sand/gravel (saturated)", 10.0, 200.0, "#E9C46A"),
            RockEntry("Granodiorite (fresh)", 1000.0, 5000.0, "#8D99AE"),
        ]
    )
    legend = GeologyLegend.from_rock_database(
        db, document_id="pcgl:test", title="Test legend"
    )
    return write_legend(legend, tmp_path / "legend.pcgl.json")


def _write_multiline_pcsf(tmp_path) -> Path:
    """Two lines, each with named stations -- MapView.from_pcsf() needs a
    station table to build 2-D curtains (unlike the point-cloud approach
    this window used before)."""
    from pycsamt.format import write_pcsf
    from pycsamt.format.multiline import build_multiline_pcsf

    profiles = {
        "L1": {
            "x": np.array([0.0, 100.0]), "z": np.array([10.0, 50.0]),
            "rho": np.array([[100.0, 110.0], [50.0, 55.0]]),
            "sta_x": np.array([0.0, 100.0]), "sta_names": ["L1-1", "L1-2"],
        },
        "L2": {
            "x": np.array([0.0, 100.0]), "z": np.array([10.0, 50.0]),
            "rho": np.array([[200.0, 210.0], [90.0, 95.0]]),
            "sta_x": np.array([0.0, 100.0]), "sta_names": ["L2-1", "L2-2"],
        },
    }
    model = build_multiline_pcsf(profiles)
    return write_pcsf(model, tmp_path / "multiline.pcsf")


def test_window_constructs(qapp):
    win = Pcsf3DWindow()
    assert win._view is None
    assert not win._btn_render.isEnabled()
    assert not win._btn_export.isEnabled()
    win.close()


@_SKIP_OCCAM
def test_load_pcsf_enables_render_and_shows_summary(qapp, tmp_path):
    win = Pcsf3DWindow()
    path = _write_grid2d_pcsf(tmp_path)
    win.load_pcsf(str(path))
    assert win._view is not None
    assert win._view.lines == ("line1",)
    assert win._btn_render.isEnabled()
    assert win._btn_export.isEnabled()
    assert "line1" not in win._file_lbl.text()  # summary uses counts, not names
    assert "lines=1" in win._file_lbl.text()
    win.close()


@_SKIP_OCCAM
def test_load_pcsf_renders_a_real_plotly_figure(qapp, tmp_path):
    win = Pcsf3DWindow()
    path = _write_grid2d_pcsf(tmp_path)
    win.load_pcsf(str(path))
    fig = win._plotly_view.figure
    assert fig is not None
    assert len(fig.data) > 0
    win.close()


@_SKIP_OCCAM
@pytest.mark.parametrize("mode", ["fence", "block", "depth", "surface"])
def test_every_mode_renders_without_error(qapp, tmp_path, mode):
    win = Pcsf3DWindow()
    path = _write_grid2d_pcsf(tmp_path)
    win.load_pcsf(str(path))
    win._combo_mode.setCurrentIndex(win._combo_mode.findData(mode))
    assert win._mode() == mode
    win._on_render()
    assert "Render error" not in win._status_lbl.text()
    assert win._plotly_view.figure is not None
    win.close()


def test_multiline_pcsf_reports_two_lines(qapp, tmp_path):
    win = Pcsf3DWindow()
    path = _write_multiline_pcsf(tmp_path)
    win.load_pcsf(str(path))
    assert win._view.lines == ("L1", "L2")
    assert win._view.n_stations == 4
    assert "lines=2" in win._file_lbl.text()
    win.close()


def test_load_pcsf_raises_and_leaves_state_untouched_for_a_bad_file(
    qapp, tmp_path
):
    win = Pcsf3DWindow()
    bad_path = tmp_path / "not_a_pcsf_file.pcsf"
    bad_path.write_text("this is not an hdf5 file")
    with pytest.raises(Exception):
        win.load_pcsf(str(bad_path))
    assert win._view is None
    assert not win._btn_render.isEnabled()
    win.close()


def test_on_load_swallows_the_error_and_shows_status(qapp, tmp_path, monkeypatch):
    win = Pcsf3DWindow()
    bad_path = tmp_path / "not_a_pcsf_file.pcsf"
    bad_path.write_text("this is not an hdf5 file")
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getOpenFileName",
        lambda *a, **k: (str(bad_path), ""),
    )
    win._on_load()
    assert win._view is None
    assert not win._btn_render.isEnabled()
    assert "error" in win._status_lbl.text().lower()
    win.close()


@_SKIP_OCCAM
def test_export_html_writes_a_file(qapp, tmp_path, monkeypatch):
    win = Pcsf3DWindow()
    path = _write_grid2d_pcsf(tmp_path)
    win.load_pcsf(str(path))

    out_path = tmp_path / "scene.html"
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getSaveFileName",
        lambda *a, **k: (str(out_path), ""),
    )
    win._on_export()

    assert out_path.exists()
    assert out_path.stat().st_size > 0
    assert "Exported to" in win._status_lbl.text()
    win.close()


def test_export_disabled_before_a_file_is_loaded(qapp):
    win = Pcsf3DWindow()
    assert not win._btn_export.isEnabled()
    win.close()


# ── Geology overlay ──────────────────────────────────────────────────────


def test_geology_controls_disabled_before_a_legend_is_loaded(qapp):
    win = Pcsf3DWindow()
    assert not win._chk_show_geology.isEnabled()
    assert not win._chk_pattern_fill.isEnabled()
    assert win._geology_kwargs() == {}
    win.close()


def test_load_legend_enables_controls_and_reports_band_count(qapp, tmp_path):
    win = Pcsf3DWindow()
    path = _write_legend_pcgl(tmp_path)
    win._legend = win._read_legend(str(path))
    win._chk_show_geology.setEnabled(True)
    win._chk_pattern_fill.setEnabled(True)
    assert len(win._legend.entries) == 2
    win.close()


def test_on_load_legend_via_dialog(qapp, tmp_path, monkeypatch):
    win = Pcsf3DWindow()
    path = _write_legend_pcgl(tmp_path)
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getOpenFileName",
        lambda *a, **k: (str(path), ""),
    )
    win._on_load_legend()
    assert win._legend is not None
    assert "2 band(s)" in win._legend_lbl.text()
    assert win._chk_show_geology.isEnabled()
    assert win._chk_pattern_fill.isEnabled()
    win.close()


def test_read_legend_accepts_a_plain_csv(qapp, tmp_path):
    win = Pcsf3DWindow()
    csv_path = tmp_path / "legend.csv"
    csv_path.write_text(
        "name,rho_min,rho_max,color\n"
        "Sand,10,200,#E9C46A\n"
        "Granite,1000,5000,#8D99AE\n"
    )
    legend = win._read_legend(str(csv_path))
    assert len(legend.entries) == 2
    win.close()


def test_on_load_legend_bad_file_shows_error(qapp, tmp_path, monkeypatch):
    win = Pcsf3DWindow()
    bad_path = tmp_path / "not_a_legend.pcgl.json"
    bad_path.write_text("not json at all")
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getOpenFileName",
        lambda *a, **k: (str(bad_path), ""),
    )
    win._on_load_legend()
    assert win._legend is None
    assert not win._chk_show_geology.isEnabled()
    assert "error" in win._legend_lbl.text().lower()
    win.close()


@_SKIP_OCCAM
def test_geology_kwargs_include_bands_when_shown(qapp, tmp_path):
    win = Pcsf3DWindow()
    win._legend = win._read_legend(str(_write_legend_pcgl(tmp_path)))
    win._chk_show_geology.setEnabled(True)
    win._chk_show_geology.setChecked(True)

    kwargs = win._geology_kwargs()
    assert "geology" in kwargs
    assert "geology_fill" not in kwargs  # pattern fill unchecked by default
    assert len(kwargs["geology"]) == 2


@_SKIP_OCCAM
def test_geology_kwargs_empty_when_overlay_unchecked(qapp, tmp_path):
    win = Pcsf3DWindow()
    win._legend = win._read_legend(str(_write_legend_pcgl(tmp_path)))
    win._chk_show_geology.setEnabled(True)
    win._chk_show_geology.setChecked(False)
    assert win._geology_kwargs() == {}


@_SKIP_OCCAM
def test_render_with_geology_overlay_applied(qapp, tmp_path):
    win = Pcsf3DWindow()
    win._legend = win._read_legend(str(_write_legend_pcgl(tmp_path)))
    win._chk_show_geology.setEnabled(True)
    win._chk_show_geology.setChecked(True)

    path = _write_grid2d_pcsf(tmp_path)
    win.load_pcsf(str(path))

    assert "Render error" not in win._status_lbl.text()
    assert win._plotly_view.figure is not None
    win.close()


@_SKIP_OCCAM
def test_render_with_pattern_fill_degrades_gracefully_without_patterns(
    qapp, tmp_path
):
    """The synthetic legend's entries carry no pattern_id, so pattern-fill
    mode must still render (empty stencil dict), not raise."""
    win = Pcsf3DWindow()
    win._legend = win._read_legend(str(_write_legend_pcgl(tmp_path)))
    win._chk_show_geology.setEnabled(True)
    win._chk_show_geology.setChecked(True)
    win._chk_pattern_fill.setEnabled(True)
    win._chk_pattern_fill.setChecked(True)

    path = _write_grid2d_pcsf(tmp_path)
    win.load_pcsf(str(path))

    assert "Render error" not in win._status_lbl.text()
    win.close()


# ── Borehole overlay ─────────────────────────────────────────────────────


def test_borehole_controls_disabled_before_a_document_is_loaded(qapp):
    win = Pcsf3DWindow()
    assert not win._chk_show_boreholes.isEnabled()
    win.close()


def test_borehole_traces_empty_without_a_document(qapp):
    win = Pcsf3DWindow()
    win._view = _geo_referenced_view()
    assert win._borehole_traces() == []
    win.close()


def test_borehole_traces_empty_when_fewer_than_two_geo_stations(qapp):
    from pycsamt.map import MapView
    from pycsamt.map._core import MapData, StationRecord

    win = Pcsf3DWindow()
    win._view = MapView(
        MapData(
            sites=None,
            stations=(StationRecord(id="S0", latitude=22.0, longitude=103.0),),
        )
    )
    win._pcbh_document = _pcbh_document()
    win._chk_show_boreholes.setEnabled(True)
    win._chk_show_boreholes.setChecked(True)
    assert win._borehole_traces() == []
    win.close()


def test_borehole_traces_built_from_a_real_geo_referenced_scene(qapp):
    win = Pcsf3DWindow()
    win._view = _geo_referenced_view()
    win._pcbh_document = _pcbh_document()
    win._chk_show_boreholes.setEnabled(True)
    win._chk_show_boreholes.setChecked(True)

    traces = win._borehole_traces()

    # 2 lithology-interval tubes (Mesh3d) + 1 collar marker (Scatter3d).
    assert len(traces) == 3
    kinds = sorted(type(t).__name__ for t in traces)
    assert kinds == ["Mesh3d", "Mesh3d", "Scatter3d"]
    win.close()


def test_borehole_traces_empty_when_overlay_unchecked(qapp):
    win = Pcsf3DWindow()
    win._view = _geo_referenced_view()
    win._pcbh_document = _pcbh_document()
    win._chk_show_boreholes.setEnabled(True)
    win._chk_show_boreholes.setChecked(False)
    assert win._borehole_traces() == []
    win.close()


def test_render_merges_borehole_traces_onto_the_scene_figure(qapp):
    win = Pcsf3DWindow()
    win._view = _geo_referenced_view()
    win._btn_render.setEnabled(True)
    win._pcbh_document = _pcbh_document()
    win._chk_show_boreholes.setEnabled(True)
    win._chk_show_boreholes.setChecked(True)

    win._on_render()

    assert "Render error" not in win._status_lbl.text()
    fig = win._plotly_view.figure
    assert fig is not None
    assert len(fig.data) == 3  # the empty synthetic scene contributes none
    win.close()


def test_on_load_boreholes_via_dialog(qapp, tmp_path, monkeypatch):
    win = Pcsf3DWindow()
    path = _write_pcbh(tmp_path)
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getOpenFileName",
        lambda *a, **k: (str(path), ""),
    )
    win._on_load_boreholes()
    assert win._pcbh_document is not None
    assert "1 hole(s)" in win._boreholes_lbl.text()
    assert win._chk_show_boreholes.isEnabled()
    win.close()


def test_on_load_boreholes_bad_file_shows_error(qapp, tmp_path, monkeypatch):
    win = Pcsf3DWindow()
    bad_path = tmp_path / "not_a_pcbh_file.pcbh.json"
    bad_path.write_text("not json at all")
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getOpenFileName",
        lambda *a, **k: (str(bad_path), ""),
    )
    win._on_load_boreholes()
    assert win._pcbh_document is None
    assert not win._chk_show_boreholes.isEnabled()
    assert "error" in win._boreholes_lbl.text().lower()
    win.close()


@_SKIP_OCCAM
def test_render_with_both_geology_and_borehole_overlays(qapp, tmp_path):
    """The two overlays are independent -- both apply in the same render."""
    win = Pcsf3DWindow()
    win._legend = win._read_legend(str(_write_legend_pcgl(tmp_path)))
    win._chk_show_geology.setEnabled(True)
    win._chk_show_geology.setChecked(True)
    win._pcbh_document = _pcbh_document()
    win._chk_show_boreholes.setEnabled(True)
    win._chk_show_boreholes.setChecked(True)

    path = _write_grid2d_pcsf(tmp_path)
    win.load_pcsf(str(path))

    assert "Render error" not in win._status_lbl.text()
    win.close()


# ── Structural overlay ───────────────────────────────────────────────────


def test_structure_controls_disabled_before_a_document_is_loaded(qapp):
    win = Pcsf3DWindow()
    assert not win._chk_show_structure.isEnabled()
    win.close()


def test_structure_traces_empty_without_a_document(qapp):
    win = Pcsf3DWindow()
    win._view = _geo_referenced_view()
    assert win._structure_traces() == []
    win.close()


def test_structure_traces_empty_when_fewer_than_two_geo_stations(qapp):
    from pycsamt.map import MapView
    from pycsamt.map._core import MapData, StationRecord

    win = Pcsf3DWindow()
    win._view = MapView(
        MapData(
            sites=None,
            stations=(StationRecord(id="S0", latitude=22.0, longitude=103.0),),
        )
    )
    win._pcgs_document = _pcgs_document()
    win._chk_show_structure.setEnabled(True)
    win._chk_show_structure.setChecked(True)
    assert win._structure_traces() == []
    win.close()


def test_structure_traces_built_from_a_real_geo_referenced_scene(qapp):
    win = Pcsf3DWindow()
    win._view = _geo_referenced_view()
    win._pcgs_document = _pcgs_document()
    win._chk_show_structure.setEnabled(True)
    win._chk_show_structure.setChecked(True)

    traces = win._structure_traces()

    # 1 fault plane (Mesh3d) + planar-measurements group + linear-measurements
    # group (both Scatter3d).
    assert len(traces) == 3
    kinds = sorted(type(t).__name__ for t in traces)
    assert kinds == ["Mesh3d", "Scatter3d", "Scatter3d"]
    win.close()


def test_structure_traces_empty_when_overlay_unchecked(qapp):
    win = Pcsf3DWindow()
    win._view = _geo_referenced_view()
    win._pcgs_document = _pcgs_document()
    win._chk_show_structure.setEnabled(True)
    win._chk_show_structure.setChecked(False)
    assert win._structure_traces() == []
    win.close()


def test_render_merges_structure_traces_onto_the_scene_figure(qapp):
    win = Pcsf3DWindow()
    win._view = _geo_referenced_view()
    win._btn_render.setEnabled(True)
    win._pcgs_document = _pcgs_document()
    win._chk_show_structure.setEnabled(True)
    win._chk_show_structure.setChecked(True)

    win._on_render()

    assert "Render error" not in win._status_lbl.text()
    fig = win._plotly_view.figure
    assert fig is not None
    assert len(fig.data) == 3  # the empty synthetic scene contributes none
    win.close()


def test_on_load_structure_via_dialog(qapp, tmp_path, monkeypatch):
    win = Pcsf3DWindow()
    path = _write_pcgs(tmp_path)
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getOpenFileName",
        lambda *a, **k: (str(path), ""),
    )
    win._on_load_structure()
    assert win._pcgs_document is not None
    assert "1 fault(s), 1 planar, 1 linear" in win._structure_lbl.text()
    assert win._chk_show_structure.isEnabled()
    win.close()


def test_on_load_structure_bad_file_shows_error(qapp, tmp_path, monkeypatch):
    win = Pcsf3DWindow()
    bad_path = tmp_path / "not_a_pcgs_file.pcgs.json"
    bad_path.write_text("not json at all")
    monkeypatch.setattr(
        "PySide6.QtWidgets.QFileDialog.getOpenFileName",
        lambda *a, **k: (str(bad_path), ""),
    )
    win._on_load_structure()
    assert win._pcgs_document is None
    assert not win._chk_show_structure.isEnabled()
    assert "error" in win._structure_lbl.text().lower()
    win.close()


def test_render_with_geology_borehole_and_structure_overlays_together(qapp):
    """All three overlays are independent -- all apply in the same render."""
    win = Pcsf3DWindow()
    win._view = _geo_referenced_view()
    win._btn_render.setEnabled(True)

    win._pcbh_document = _pcbh_document()
    win._chk_show_boreholes.setEnabled(True)
    win._chk_show_boreholes.setChecked(True)

    win._pcgs_document = _pcgs_document()
    win._chk_show_structure.setEnabled(True)
    win._chk_show_structure.setChecked(True)

    win._on_render()

    assert "Render error" not in win._status_lbl.text()
    fig = win._plotly_view.figure
    # 3 borehole traces + 3 structure traces, no data traces of its own.
    assert len(fig.data) == 6
    win.close()


# ── scene controls, topography sources, Map View (2026-09-25) ────────────


@pytest.fixture
def loaded(qapp, tmp_path):
    win = Pcsf3DWindow()
    win.load_pcsf(str(_write_multiline_pcsf(tmp_path)))
    yield win
    win.close()


def test_scene_is_white_with_exaggeration(loaded):
    L = loaded._plotly_view.figure.layout
    assert L.paper_bgcolor == "#ffffff" and L.scene.bgcolor == "#ffffff"
    assert L.scene.aspectmode == "manual"
    assert "vertical ×" in loaded._status_lbl.text()


def test_colorbar_legend_and_background_toggles(loaded):
    loaded._chk_colorbar.setChecked(False)
    loaded._chk_legend.setChecked(False)
    loaded._combo_bg.setCurrentText("Dark")
    fig = loaded._plotly_view.figure
    assert all(not tr.showscale for tr in fig.data
               if "showscale" in tr and tr.showscale is not None)
    assert fig.layout.showlegend is False
    assert fig.layout.paper_bgcolor == "#111827"


def test_manual_exaggeration_scales_depth(loaded):
    loaded._spin_ve.setValue(1.0)
    loaded._on_render()
    z1 = loaded._plotly_view.figure.layout.scene.aspectratio.z
    loaded._spin_ve.setValue(2.0)
    loaded._on_render()
    z2 = loaded._plotly_view.figure.layout.scene.aspectratio.z
    assert z2 == pytest.approx(min(2 * z1, 3.0))


def test_topography_from_a_name_elevation_table(loaded, tmp_path):
    names = [s.id for s in loaded._view.data.stations]
    p = tmp_path / "topo.csv"
    p.write_text("station,elevation\n" + "\n".join(
        f"{n},{500 + i}" for i, n in enumerate(names)))
    assert loaded.load_topography(str(p)) == len(names)
    assert loaded._topo_source() == "topo_file"
    assert f"{len(names)}/{len(names)}" in loaded._topo_lbl.text()


def test_topography_from_loaded_stations_needs_a_survey(loaded):
    loaded._combo_topo.setCurrentIndex(
        loaded._combo_topo.findData("stations"))
    assert "No survey loaded" in loaded._topo_lbl.text()


def test_flat_topography(loaded):
    loaded._combo_topo.setCurrentIndex(loaded._combo_topo.findData("flat"))
    assert "Flat" in loaded._topo_lbl.text()


def test_spin_toggles_the_page_flag(loaded, monkeypatch):
    calls = []
    monkeypatch.setattr(loaded._plotly_view, "run_js", calls.append)
    loaded._btn_spin.setChecked(True)
    loaded._btn_spin.setChecked(False)
    assert calls == ["window.__pycsamtSpin = true;",
                     "window.__pycsamtSpin = false;"]


def test_open_in_map_view_launches_with_the_file(loaded, monkeypatch):
    import pycsamt.app.desktop.mapview_bridge as mb

    seen = {}

    def fake(pcsf, **kw):
        seen["pcsf"] = pcsf
        return mb.MapViewLaunch(url="http://127.0.0.1:8771")

    monkeypatch.setattr(mb, "launch_mapview", fake)
    assert loaded._btn_mapview.isEnabled()
    loaded._on_open_mapview()
    assert seen["pcsf"].endswith(".pcsf")
    assert "8771" in loaded._status_lbl.text()


def test_mapview_command_frozen_and_source(monkeypatch):
    import sys

    import pycsamt.app.desktop.mapview_bridge as mb

    cmd = mb.mapview_command("m.pcsf", "127.0.0.1", 8775)
    assert cmd[1:3] == ["-m", "pycsamt.app.mapview"]
    assert cmd[-2:] == ["--pcsf", "m.pcsf"] and "8775" in cmd
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    assert mb.mapview_command(None, "h", 1)[1] == mb.MAPVIEW_SERVER_FLAG


# ── v2.6 toolbar, topography and Map View hand-off ─────────────────────────

_DEMO = Path(__file__).resolve().parents[4] / "examples"
_TOPO_PCSF = _DEMO / "pcsf_conversion_demo" / "output" / "occam2d_with_topo.pcsf"
_RASTER_PCSM = (_DEMO / "pcsm_conversion_demo" / "output"
                / "modem3d_with_raster_topo.pcsm")


def test_toolbar_offers_jet_r_and_map_view_colours(loaded):
    from pycsamt.app.desktop.controllers.pcsf_scene import CMAPS

    items = [loaded._combo_cmap.itemText(i)
             for i in range(loaded._combo_cmap.count())]
    assert "jet_r" in items and items == list(CMAPS)
    loaded._combo_cmap.setCurrentText("jet_r")
    assert "jet_r" in loaded._status_lbl.text()


def test_mode_dropdown_drives_the_scene(loaded):
    labels = [loaded._combo_mode.itemText(i)
              for i in range(loaded._combo_mode.count())]
    assert labels == ["Fence", "Block", "Depth slices", "Iso-surface"]
    loaded._combo_mode.setCurrentIndex(loaded._combo_mode.findData("depth"))
    assert loaded._status_lbl.text().startswith("depth")


def test_toolbar_fits_a_1280_px_screen(loaded):
    # 1920 x 1080 at 150 % = 1280 logical px; the scene toolbar must fit
    # beside the side panel (and fully when the panel is hidden)
    from PySide6.QtGui import QFontDatabase

    bar = loaded.findChild(QWidget, "SceneToolbar")
    # one row of controls (the old two-row toolbar was ~1500 px wide)
    assert bar.minimumSizeHint().height() <= 2 *         loaded._combo_cmap.sizeHint().height()
    if QFontDatabase.families():  # pixel widths are only real with fonts
        assert bar.minimumSizeHint().width() <= 1280 - 250 - 20
    loaded._btn_panel.setChecked(False)
    assert not loaded._splitter.widget(0).isVisible()
    loaded._btn_panel.setChecked(True)


def test_layers_menu_toggles(loaded):
    acts = loaded._btn_layers.menu().actions()
    names = [a.text().split("	")[0] for a in acts if a.isCheckable()]
    assert names == ["Topography", "Terrain line", "Colour bar", "Legend",
                     "Station markers", "Station labels"]
    loaded._chk_colorbar.setChecked(False)
    fig = loaded._plotly_view.figure
    assert all(not tr.showscale for tr in fig.data
               if "showscale" in tr and tr.showscale is not None)


def test_depth_preset_and_custom(loaded):
    i = loaded._combo_depth.findData(1000.0)
    loaded._combo_depth.setCurrentIndex(i)
    assert loaded._depth_max() == 1000.0
    assert "to 1,000 m" in loaded._status_lbl.text()
    loaded._combo_depth.setCurrentIndex(loaded._combo_depth.findData(-1.0))
    assert loaded._spin_depth.isVisibleTo(loaded)


def test_topography_and_terrain_toggles(loaded):
    loaded._chk_topo.setChecked(False)
    assert loaded._view is not None
    loaded._chk_terrain.setChecked(False)
    assert "Render error" not in loaded._status_lbl.text()


def test_camera_presets_run_in_the_page(loaded, monkeypatch):
    calls = []
    monkeypatch.setattr(loaded._plotly_view, "run_js", calls.append)
    loaded.set_camera("top")
    assert "scene.camera" in calls[-1] and '"z": 2.5' in calls[-1]
    assert loaded._camera["eye"]["z"] == 2.5
    # the camera is kept on the next render and handed to Map View
    loaded._on_render()
    assert loaded._plotly_view.figure.layout.scene.camera.eye.z == 2.5
    assert loaded.mapview_state()["camera"]["eye"]["z"] == 2.5


def test_step_exaggeration(loaded):
    loaded._spin_ve.setValue(2.0)
    loaded._step_ve(+1)
    assert loaded._spin_ve.value() == 3.0
    loaded._step_ve(-1)
    assert loaded._spin_ve.value() == 2.0


def test_shortcuts_installed(loaded):
    from PySide6.QtGui import QShortcut

    keys = {s.key().toString() for s in loaded.findChildren(QShortcut)}
    assert {"1", "4", "Space", "R", "T", "C", "L", "S", "F5",
            "Ctrl+M", "Ctrl+P"} <= keys


@pytest.mark.skipif(not _TOPO_PCSF.is_file(), reason="demo PCSF missing")
def test_auto_topography_from_the_file_fills_gaps(qapp):
    win = Pcsf3DWindow()
    win.load_pcsf(str(_TOPO_PCSF))  # 2 of 47 stations have an elevation
    assert win._topo_source() == "auto"
    assert "station table" in win._file_lbl.text()
    assert "45 filled along lines" in win._topo_lbl.text()
    assert len(win._elev_used) == 47
    win._chk_fill_gaps.setChecked(False)
    assert "(2/47 stations)" in win._topo_lbl.text()
    win.close()


@pytest.mark.skipif(not _RASTER_PCSM.is_file(), reason="demo PCSM missing")
def test_file_elevations_sample_a_topography_raster():
    from pycsamt.app.desktop.controllers.pcsf_scene import file_elevations

    elev, how = file_elevations(_RASTER_PCSM)
    assert "station table" in how
    assert len(elev) >= 112
    assert all(abs(e) < 1e4 for e in elev.values())


def test_fill_line_gaps_interpolates_between_known_stations():
    from pycsamt.app.desktop.controllers.pcsf_scene import fill_line_gaps
    from pycsamt.map._core import StationRecord

    st = [StationRecord(f"S{i}", 0.0, float(i), line="L1")
          for i in range(5)]
    out, n = fill_line_gaps(st, {"S0": 100.0, "S4": 140.0})
    assert n == 3 and out["S2"] == pytest.approx(120.0)
    out, n = fill_line_gaps(st, {"S0": 100.0})  # one point: nothing to do
    assert n == 0


def test_mapview_state_carries_the_scene(loaded):
    loaded._combo_cmap.setCurrentText("jet_r")
    loaded._combo_mode.setCurrentIndex(loaded._combo_mode.findData("block"))
    loaded._btn_spin.setChecked(True)
    st = loaded.mapview_state()
    c = st["controls"]
    assert st["view"] == "map3d" and st["spin"] is True
    assert c["mode3d"] == "block" and c["cmap"] == "jet_r"
    assert c["show_stations"] is True and "elevations" in st
    import json

    json.dumps(st)  # JSON-safe: goes to Map View as --state


def test_open_in_map_view_hands_over_the_scene(loaded, monkeypatch):
    import pycsamt.app.desktop.mapview_bridge as mb

    seen = {}

    def fake(pcsf, **kw):
        seen.update(kw, pcsf=pcsf)
        return mb.MapViewLaunch(url="http://127.0.0.1:8772")

    monkeypatch.setattr(mb, "launch_mapview", fake)
    loaded._combo_depth.setCurrentIndex(loaded._combo_depth.findData(500.0))
    loaded._on_open_mapview()
    assert seen["state"]["controls"]["depth_hi"] == 500.0
    assert "with this scene" in loaded._status_lbl.text()


def test_bridge_writes_state_file(tmp_path):
    import json

    import pycsamt.app.desktop.mapview_bridge as mb

    path = mb.write_state({"view": "map3d"})
    assert json.loads(Path(path).read_text())["view"] == "map3d"
    cmd = mb.mapview_command("m.pcsf", "127.0.0.1", 8775, path)
    assert cmd[-2:] == ["--state", path]


def test_mapview_seed_state_outputs():
    from dash import no_update

    from pycsamt.app.mapview.callbacks.chrome import (
        _SEED_WIDGETS,
        seed_state_outputs,
    )

    assert all(o is no_update for o in seed_state_outputs(None))
    out = seed_state_outputs({"view": "map3d", "spin": True,
                              "controls": {"mode3d": "depth",
                                           "cmap": "jet_r"},
                              "camera": {"eye": {"x": 1, "y": 1, "z": 1}}})
    keys = [k for _w, k in _SEED_WIDGETS]
    assert out[keys.index("mode3d")] == "depth"
    assert out[keys.index("cmap")] == "jet_r"
    assert out[keys.index("depth_hi")] is no_update
    n = len(_SEED_WIDGETS)
    assert out[n + 3] is True  # spin
    assert out[n + 5] == {"map3d": {"camera": {"eye": {"x": 1, "y": 1,
                                                       "z": 1}}}}
    assert out[n + 6] == 1  # clicks the 3-D rail


def test_mapview_main_applies_state_elevations(tmp_path, monkeypatch):
    import json
    import sys

    import pycsamt.app.mapview.__main__ as mvmain
    import pycsamt.app.mapview.app as mvapp

    if not _TOPO_PCSF.is_file():
        pytest.skip("demo PCSF missing")
    st = tmp_path / "s.json"
    st.write_text(json.dumps({"view": "map3d",
                              "elevations": {"S05": 999.0}}))
    got = {}
    monkeypatch.setattr(mvapp, "launch", lambda **kw: got.update(kw))
    monkeypatch.setattr(sys, "argv", ["x", "--pcsf", str(_TOPO_PCSF),
                                      "--state", str(st), "--no-browser"])
    assert mvmain.main() == 0
    el = {s.id: s.elevation for s in got["view"].data.stations}
    assert el["S05"] == 999.0 and got["state"]["view"] == "map3d"


# ── render controls panel (stations, labels, resistivity, view) ────────────

_BLN_PCSF = _DEMO / "pcsf_conversion_demo" / "output" / "occam2d_with_bln_topo.pcsf"


@pytest.fixture
def real(qapp):
    if not _BLN_PCSF.is_file():
        pytest.skip("demo PCSF missing")
    win = Pcsf3DWindow()
    win.load_pcsf(str(_BLN_PCSF))
    yield win
    win.close()


def test_controls_panel_hidden_by_default(qapp):
    win = Pcsf3DWindow()
    win.show()
    assert not win._btn_panel.isChecked()
    assert not win._splitter.widget(0).isVisible()
    win._btn_panel.setChecked(True)
    assert win._splitter.widget(0).isVisible()
    assert {"Stations & labels", "Depth & units", "Resistivity & colours",
            "View & geometry", "Topography"} <= set(win._sections)
    win.close()


def test_custom_depth_opens_the_panel(loaded):
    loaded._combo_depth.setCurrentIndex(loaded._combo_depth.findData(-1.0))
    assert loaded._btn_panel.isChecked()
    assert loaded._sections["Depth & units"].isChecked()


def test_station_marker_options_reach_the_scene(real):
    real._combo_sta_symbol.setCurrentIndex(
        real._combo_sta_symbol.findData("triangle-down"))
    real.set_station_color("#ff0000")
    real._rerender()
    marks = [t for t in real._plotly_view.figure.data
             if t.type == "mesh3d" and t.name == "stations"]
    assert marks and marks[0].color == "#ff0000"


def test_label_rotation_and_density(real):
    real._chk_station_labels.setChecked(True)
    real._spin_label_angle.setValue(45)
    real._combo_label_density.setCurrentIndex(
        real._combo_label_density.findData(0.5))
    real._rerender()
    ann = real._plotly_view.figure.layout.scene.annotations
    assert ann and all(a.textangle == 45 for a in ann)
    assert len(ann) < real._view.n_stations  # 50 % labelled


def test_resistivity_colour_range(real):
    real._spin_vmin.setValue(10)
    real._spin_vmax.setValue(1000)
    real._rerender()
    surf = [t for t in real._plotly_view.figure.data if t.type == "surface"]
    assert surf[0].cmin == pytest.approx(1.0)  # log10(10)
    assert surf[0].cmax == pytest.approx(3.0)


def test_proportions_modes(real):
    real._rerender()
    assert real._plotly_view.figure.layout.scene.aspectmode == "manual"
    real._combo_aspect.setCurrentIndex(real._combo_aspect.findData("data"))
    real._rerender()
    assert real._plotly_view.figure.layout.scene.aspectmode == "data"


def test_reset_render_settings(real):
    real._spin_label_angle.setValue(60)
    real._spin_vmin.setValue(5)
    real.set_station_color("#00ff00")
    real.reset_render_settings()
    opts = real.render_options()
    assert opts["station_label_angle"] == 0.0
    assert opts["value_range"] is None
    assert opts["station_color"] == "#111111"
    assert opts["station_symbol"] == "triangle-down"


def test_render_options_to_map_view_controls(real):
    from pycsamt.app.desktop.controllers.pcsf_scene import mapview_controls

    real._combo_label_density.setCurrentIndex(
        real._combo_label_density.findData(0.25))
    real._combo_scale.setCurrentIndex(real._combo_scale.findData("linear"))
    real._spin_sta_max.setValue(10)
    real._edit_label_names.setText("S01, S05")
    c = mapview_controls(real.render_options())
    assert c["station_label_density"] == "0.25"
    assert c["scale"] == "linear"
    assert c["station_max"] == 10
    assert c["station_label_names"] == "S01, S05"
    assert c["aspect"] == "data" and c["vertical_exaggeration"] == 0.0
    assert c["section_res"] == "100" and c["volume_smoothing"] == "0"


def test_vertical_exaggeration_keeps_a_single_line_box():
    import plotly.graph_objects as go

    from pycsamt.map.volume import apply_vertical_exaggeration

    fig = go.Figure(go.Scatter3d(x=[0, 10000], y=[0, 0], z=[0, -2000]))
    ve = apply_vertical_exaggeration(fig, 0)
    r = fig.layout.scene.aspectratio
    assert r.y == pytest.approx(0.2) and r.x == pytest.approx(1.0)
    assert ve == pytest.approx(2.5)  # depth half the length


def test_map_view_has_vertical_exaggeration_control():
    from pycsamt.app.mapview._ids import IDs
    from pycsamt.app.mapview.callbacks.chrome import _SEED_WIDGETS

    assert (IDs.CTL_VE, "vertical_exaggeration") in _SEED_WIDGETS
    assert (IDs.CTL_STA_SYMBOL, "station_symbol") in _SEED_WIDGETS
    from pycsamt.app.mapview.app import create_app

    app = create_app()
    key = next(k for k in app.callback_map
               if k.startswith(IDs.STORE_CONTROLS)
               or f".{IDs.STORE_CONTROLS}." in k)
    inputs = [i["id"] for i in app.callback_map[key]["inputs"]]
    assert IDs.CTL_VE in inputs


def test_scene_choices_live_in_the_controls_panel(loaded):
    panel = loaded._splitter.widget(0)
    for w in (loaded._combo_mode, loaded._combo_depth, loaded._combo_cmap,
              loaded._spin_ve):
        assert panel.isAncestorOf(w)
    assert loaded._sections["Scene"].isChecked()  # open by default


def test_top_bar_is_icon_shortcuts_only(loaded):
    from PySide6.QtWidgets import QComboBox, QDoubleSpinBox

    bar = loaded.findChild(QWidget, "SceneToolbar")
    assert not bar.findChildren(QComboBox)
    assert not bar.findChildren(QDoubleSpinBox)
    for b in (loaded._btn_panel, loaded._btn_open, loaded._btn_layers,
              loaded._btn_camera, loaded._btn_spin, loaded._btn_out,
              loaded._btn_mapview):
        assert b.toolTip()


def test_hard_render_button_on_the_scene(loaded, monkeypatch):
    assert loaded._btn_render.parent() is loaded._plotly_view
    assert loaded._btn_render.isEnabled()
    calls = []
    monkeypatch.setattr(loaded._plotly_view, "reload_next",
                        lambda: calls.append(1))
    before = loaded._view
    loaded._btn_render.click()
    assert calls == [1]
    assert loaded._view is not before  # the file was read again
    assert "Render error" not in loaded._status_lbl.text()
