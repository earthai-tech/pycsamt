# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for MapPanel.load_pcsf() / MapViewerWindow's "Load PCSF stations…".

A PCSF/PCSM file is an inversion model: since 2026-09-25 it drives the
Depth map (true model resistivity at a depth, from the per-station
inversion sections of ``pycsamt.map.MapView.from_pcsf``) and the
Resistivity map (log-mean model resistivity from the surface to that
depth).  Stations without lon/lat are geo-referenced by name against the
survey already on the map, else plotted in the model's local x/y.
"""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.panels.map_panel import MapPanel
from pycsamt.app.desktop.windows.map_window import MapViewerWindow


def _write_geo_multiline_pcsf(tmp_path) -> Path:
    """A real, round-tripped multiline PCSF with two geo-referenced
    stations -- the shape MapView.from_pcsf() needs to build a station
    table (see fix_station_keys_accept_numpy_arrays regression: sta_x/
    sta_names/sta_lat/sta_lon/sta_elev all accept numpy arrays)."""
    from pycsamt.format import write_pcsf
    from pycsamt.format.multiline import build_multiline_pcsf

    profiles = {
        "L1": {
            "x": np.array([0.0, 100.0]),
            "z": np.array([10.0, 50.0]),
            "rho": np.array([[100.0, 110.0], [50.0, 55.0]]),
            "sta_x": np.array([0.0, 100.0]),
            "sta_names": ["L1-1", "L1-2"],
            "sta_lat": np.array([22.000, 22.001]),
            "sta_lon": np.array([103.000, 103.001]),
            "sta_elev": np.array([500.0, 505.0]),
        },
    }
    model = build_multiline_pcsf(profiles)
    return write_pcsf(model, tmp_path / "geo_multiline.pcsf")


def _write_grid2d_pcsf(tmp_path) -> Path:
    """A real Occam2D-derived PCSF -- profile-only, no real lat/lon
    unless geo-referenced against an EDI survey (has_geo=False)."""
    from pycsamt.format import write_pcsf
    from pycsamt.format.adapters.occam2d import occam2d_to_pcsf
    from pycsamt.models.occam2d.results import InversionResult

    _ROOT = Path(__file__).parents[4]
    occam_dir = _ROOT / "data" / "occam2D"
    if not occam_dir.exists():
        pytest.skip(f"bundled occam2D data not found: {occam_dir}")
    result = InversionResult(workdir=occam_dir)
    model = occam2d_to_pcsf(result)
    return write_pcsf(model, tmp_path / "grid2d.pcsf")


# ── MapPanel.load_pcsf() ─────────────────────────────────────────────────


def test_load_pcsf_populates_dataframe_with_real_geo_stations(qapp, tmp_path):
    panel = MapPanel()
    path = _write_geo_multiline_pcsf(tmp_path)
    panel.load_pcsf(str(path))

    assert list(panel._df["ID"]) == ["L1-1", "L1-2"]
    np.testing.assert_allclose(panel._df["Latitude"], [22.000, 22.001])
    np.testing.assert_allclose(panel._df["Longitude"], [103.000, 103.001])
    np.testing.assert_allclose(panel._df["Elevation"], [500.0, 505.0])
    panel.close()


def test_load_pcsf_clears_stale_sites(qapp, tmp_path):
    """A prior EDI-sourced _sites must not linger against unrelated PCSF
    station IDs (Depth/Resistivity would silently use the wrong survey)."""
    panel = MapPanel()
    panel._sites = object()  # stand-in for a previously loaded EDI Sites
    path = _write_geo_multiline_pcsf(tmp_path)
    panel.load_pcsf(str(path))
    assert panel._sites is None
    panel.close()


def test_load_pcsf_without_geo_uses_local_model_coordinates(
    qapp, tmp_path
):
    """An Occam2D result has no lat/lon: its stations are placed on the
    model's local x (metres) instead of being dropped (they used to
    vanish, so the map showed nothing)."""
    panel = MapPanel()
    path = _write_grid2d_pcsf(tmp_path)
    panel.load_pcsf(str(path))
    assert panel._pcsf_local and len(panel._df) == 47
    assert panel._df["Longitude"].is_monotonic_increasing  # profile x
    panel.close()


def test_load_pcsf_is_georeferenced_against_the_survey_on_the_map(
    qapp, tmp_path
):
    import pandas as pd

    panel = MapPanel()
    names = [f"S{i:02d}" for i in range(47)]
    panel.set_dataframe(pd.DataFrame({
        "ID": names, "Latitude": np.linspace(26.0, 26.01, 47),
        "Longitude": np.linspace(110.0, 110.02, 47)}))
    panel.load_pcsf(str(_write_grid2d_pcsf(tmp_path)))
    assert not panel._pcsf_local
    np.testing.assert_allclose(panel._df["Latitude"].iloc[[0, -1]],
                               [26.0, 26.01])
    panel.close()


def test_model_depth_and_mean_resistivity_maps(qapp, tmp_path):
    from pycsamt.map._core import resistivity_at_depth

    panel = MapPanel()
    panel.load_pcsf(str(_write_geo_multiline_pcsf(tmp_path)))
    assert panel.has_model
    assert panel.model_depth_range() == (10.0, 50.0)
    at30 = panel._model_rho_at_depth(30.0)
    assert at30 == resistivity_at_depth(panel._pcsf_view.data, 30.0)
    assert at30["L1-1"] == pytest.approx(75.0)  # halfway 100 -> 50
    mean = panel._model_rho_average(50.0)
    assert mean["L1-1"] == pytest.approx(np.sqrt(100.0 * 50.0))
    for mt in ("depth", "resistivity"):
        panel._map_type = mt
        panel._target_depth_m = 30.0
        panel._draw_map()
        assert panel._scatter is not None
        assert "model" in panel._canvas.figure.axes[0].get_title().lower() \
            or "Model" in panel._canvas.figure.axes[0].get_title()
    panel._map_type = "depth"
    panel._target_depth_m = 5000.0
    panel._draw_map()
    assert "outside the model" in panel._canvas.figure.axes[0].get_title()
    panel.close()


def test_loading_survey_data_leaves_model_mode(qapp, tmp_path):
    panel = MapPanel()
    panel.load_pcsf(str(_write_geo_multiline_pcsf(tmp_path)))
    panel.set_sites([])
    assert not panel.has_model and not panel._pcsf_local
    panel.close()


def test_load_pcsf_fetch_elevation_defaults_off(qapp, tmp_path):
    """fetch_elevation defaults False -- must not attempt a network call."""
    panel = MapPanel()
    path = _write_geo_multiline_pcsf(tmp_path)
    with mock.patch("pycsamt.map.MapView.from_pcsf") as m:
        m.return_value = mock.Mock(
            data=mock.Mock(stations=[])
        )
        panel.load_pcsf(str(path))
        assert m.call_args.kwargs["fetch_elevation"] is False
    panel.close()


def test_load_pcsf_bad_file_raises(qapp, tmp_path):
    panel = MapPanel()
    bad_path = tmp_path / "not_a_pcsf_file.pcsf"
    bad_path.write_text("this is not an hdf5 file")
    with pytest.raises(Exception):
        panel.load_pcsf(str(bad_path))
    panel.close()


# ── MapViewerWindow "Load PCSF stations…" ────────────────────────────────


def test_window_load_pcsf_via_dialog(qapp, tmp_path):
    win = MapViewerWindow()
    path = _write_geo_multiline_pcsf(tmp_path)
    with mock.patch(
        "PySide6.QtWidgets.QFileDialog.getOpenFileName",
        return_value=(str(path), ""),
    ):
        win._on_load_pcsf()

    assert "2 station(s)" in win._lbl_pcsf.text()
    assert "model depth 10–50 m" in win._lbl_pcsf.text()
    assert win._map_panel._df.shape[0] == 2
    # A model opens on the Depth map, depth box clamped to the model,
    # no frequency box (a model has depths, not frequencies).
    assert win._combo_type.currentText() == "Depth"
    assert (win._spin_depth.minimum(), win._spin_depth.maximum()) == (
        10.0, 50.0)
    assert win._grp_freq.isHidden()
    win._combo_type.setCurrentText("Resistivity")
    assert not win._grp_depth.isHidden()
    assert win._grp_depth.title() == "Average from the surface to"
    win.close()


def test_window_load_pcsf_dialog_cancelled_is_a_noop(qapp):
    win = MapViewerWindow()
    with mock.patch(
        "PySide6.QtWidgets.QFileDialog.getOpenFileName",
        return_value=("", ""),
    ):
        win._on_load_pcsf()
    assert win._lbl_pcsf.text() == ""
    win.close()


def test_window_load_pcsf_bad_file_shows_error(qapp, tmp_path):
    win = MapViewerWindow()
    bad_path = tmp_path / "not_a_pcsf_file.pcsf"
    bad_path.write_text("this is not an hdf5 file")
    with mock.patch(
        "PySide6.QtWidgets.QFileDialog.getOpenFileName",
        return_value=(str(bad_path), ""),
    ):
        win._on_load_pcsf()
    assert "error" in win._lbl_pcsf.text().lower()
    win.close()


def test_window_load_pcsf_fetch_elevation_checkbox_forwarded(qapp, tmp_path):
    win = MapViewerWindow()
    win._chk_pcsf_fetch_elev.setChecked(True)
    path = _write_geo_multiline_pcsf(tmp_path)
    with mock.patch.object(
        win._map_panel, "load_pcsf", wraps=win._map_panel.load_pcsf
    ) as m:
        with mock.patch(
            "PySide6.QtWidgets.QFileDialog.getOpenFileName",
            return_value=(str(path), ""),
        ):
            win._on_load_pcsf()
        assert m.call_args.kwargs["fetch_elevation"] is True
    win.close()


# ── 3-D model slice + library map types (2026-09-25) ─────────────────────


def _bh_pcsf(tmp_path) -> Path:
    bh = Path(__file__).parents[4] / "data" / "MT" / "broken-hill" / \
        "final-models"
    if not bh.is_dir():
        pytest.skip("Broken Hill ModEM result not present")
    from pycsamt.format import convert_engine as ce

    out = tmp_path / "bh.pcsf"
    ce.write_model(ce.build_model(ce.detect(bh, None), ce.ConvertOptions()),
                   out, "pcsf", False)
    return out


def test_grid3d_depth_map_draws_the_model_slice(qapp, tmp_path):
    from matplotlib.collections import QuadMesh

    panel = MapPanel()
    panel.load_pcsf(str(_bh_pcsf(tmp_path)))
    assert panel._pcsf_kind == "grid3d" and panel._pcsf_model is not None
    panel._map_type = "depth"
    panel._target_depth_m = 500.0
    panel._draw_map()
    ax = panel._canvas.figure.axes[0]
    assert any(isinstance(c, QuadMesh) for c in ax.collections)
    # station markers share the slice's colour limits
    mesh = next(c for c in ax.collections if isinstance(c, QuadMesh))
    assert panel._scatter.get_clim() == mesh.get_clim()
    panel.set_show_model_slice(False)
    ax = panel._canvas.figure.axes[0]
    assert not any(isinstance(c, QuadMesh) for c in ax.collections)
    panel.close()


@pytest.fixture(scope="module")
def kap_sites():
    root = Path(__file__).parents[4] / "data" / "MT" / "kap03lmt_edis"
    if not root.is_dir():
        pytest.skip("kap03 EDIs not present")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(root))


@pytest.mark.parametrize("mt", ["phase tensor", "induction arrows",
                                "strike", "dimensionality", "confidence"])
def test_emtools_map_types_draw_from_sites(qapp, kap_sites, mt):
    panel = MapPanel()
    panel.set_sites(kap_sites)
    panel.redraw(map_type=mt, target_freq_hz=0.01)
    ax = panel._canvas.figure.axes[0]
    assert ax.get_title()
    assert not any("could not be drawn" in t.get_text() for t in ax.texts)
    assert ax.collections or ax.patches or ax.lines
    panel.close()


def test_emtools_map_explains_missing_data_for_a_model(qapp, tmp_path):
    panel = MapPanel()
    panel.load_pcsf(str(_write_geo_multiline_pcsf(tmp_path)))
    panel.redraw(map_type="phase tensor")
    texts = [t.get_text() for t in panel._canvas.figure.axes[0].texts]
    assert any("no transfer" in t for t in texts)
    panel.close()


def test_period_snaps_to_the_survey(qapp, kap_sites):
    panel = MapPanel()
    panel.set_sites(kap_sites)
    panel._target_freq_hz = 0.012
    f = panel._nearest_survey_freq()
    assert f in panel._collect_frequencies()
    panel.close()
