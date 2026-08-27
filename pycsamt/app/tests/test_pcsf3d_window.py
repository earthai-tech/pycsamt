# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for Pcsf3DWindow (pycsamt.app.desktop.windows.pcsf3d_window).

Phase 7 of the PCSF format plan: the desktop app's first 3-D/volume
panel. Loads a real Occam2D-derived .pcsf file (built from the real
bundled data/occam2D dataset) so the rendered scatter is genuine
inverted resistivity, not synthetic placeholder data.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.windows.pcsf3d_window import Pcsf3DWindow

_ROOT = Path(__file__).parents[3]
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


def _write_multiline_pcsf(tmp_path) -> Path:
    from pycsamt.format import write_pcsf
    from pycsamt.format.multiline import build_multiline_pcsf

    profiles = {
        "L1": {"x": np.array([0.0, 100.0]), "z": np.array([10.0, 50.0]),
               "rho": np.array([[100.0, 110.0], [50.0, 55.0]])},
        "L2": {"x": np.array([0.0, 100.0]), "z": np.array([10.0, 50.0]),
               "rho": np.array([[200.0, 210.0], [90.0, 95.0]])},
    }
    model = build_multiline_pcsf(profiles)
    return write_pcsf(model, tmp_path / "multiline.pcsf")


def test_window_constructs(qapp):
    win = Pcsf3DWindow()
    assert win._model is None
    assert not win._btn_render.isEnabled()
    win.close()


@_SKIP_OCCAM
def test_load_pcsf_enables_render_and_shows_summary(qapp, tmp_path):
    win = Pcsf3DWindow()
    path = _write_grid2d_pcsf(tmp_path)
    win.load_pcsf(str(path))
    assert win._model is not None
    assert win._model.kind == "grid2d"
    assert win._btn_render.isEnabled()
    assert "grid2d" in win._file_lbl.text()
    win.close()


@_SKIP_OCCAM
def test_load_pcsf_renders_a_real_point_cloud(qapp, tmp_path):
    win = Pcsf3DWindow()
    path = _write_grid2d_pcsf(tmp_path)
    win.load_pcsf(str(path))
    assert win._cloud is not None
    assert win._cloud.x.size > 0
    assert np.all(np.isfinite(win._cloud.value))
    # A real Axes3D was actually built on the shared canvas figure.
    assert win._canvas.figure.axes
    assert win._canvas.figure.axes[0].name == "3d"
    win.close()


def test_multiline_pcsf_renders_two_offset_lines(qapp, tmp_path):
    win = Pcsf3DWindow()
    path = _write_multiline_pcsf(tmp_path)
    win.load_pcsf(str(path))
    assert win._model.kind == "multiline"
    assert "lines=2" in win._file_lbl.text()
    assert set(np.unique(win._cloud.y).tolist()) == {0.0, 1000.0}
    win.close()


def test_max_points_spinbox_limits_render(qapp, tmp_path):
    from pycsamt.format import write_pcsf
    from pycsamt.format.multiline import build_multiline_pcsf

    # 2 lines x 50x30 points = 3000 total, well above the spinbox's own
    # UI-realistic minimum (1000), so the cap is actually exercised.
    x = np.linspace(0, 1000, 50)
    z = np.linspace(1, 500, 30)
    rho = np.full((30, 50), 100.0)
    profiles = {"L1": {"x": x, "z": z, "rho": rho}, "L2": {"x": x, "z": z, "rho": rho}}
    model = build_multiline_pcsf(profiles, cache_derived_volume=False)
    path = write_pcsf(model, tmp_path / "big_multiline.pcsf")

    win = Pcsf3DWindow()
    win._spin_max_points.setValue(1000)
    win.load_pcsf(str(path))
    assert win._cloud.x.size == 1000
    win.close()


def test_load_pcsf_raises_and_leaves_state_untouched_for_a_bad_file(
    qapp, tmp_path
):
    win = Pcsf3DWindow()
    bad_path = tmp_path / "not_a_pcsf_file.pcsf"
    bad_path.write_text("this is not an hdf5 file")
    with pytest.raises(Exception):
        win.load_pcsf(str(bad_path))
    assert win._model is None
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
    assert win._model is None
    assert not win._btn_render.isEnabled()
    assert "error" in win._status_lbl.text().lower()
    win.close()
