# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the Airborne EM studio (controllers/airborne_studio +
windows/airborne_window), on the bundled synthetic EMTF-XML surveys.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers import airborne_studio as st
from pycsamt.app.desktop.controllers.correction_views import (
    figure_blank_reason,
)
from pycsamt.app.desktop.controllers.inversion_engines import defaults

_ROOT = Path(__file__).resolve().parents[4]
DATA = {
    "ztem": _ROOT / "data" / "ZTEM" / "forrestania_wa",
    "afmag_original": _ROOT / "data" / "AFMAG" / "abitibi_on",
    "afmag_airmt": _ROOT / "data" / "AFMAG" / "yulong_belt_cn",
    "mobilemt": _ROOT / "data" / "mobileMT" / "timiskaming_kimberlite_on",
}
_HAVE = {k: p.exists() and any(p.rglob("*.xml")) for k, p in DATA.items()}


def _need(tech):
    return pytest.mark.skipif(not _HAVE[tech], reason=f"{tech} data missing")


@pytest.fixture(scope="module")
def surveys():
    from pycsamt.airborne.site import ensure_asites

    return {k: ensure_asites(str(p), recursive=True, strict=False, verbose=0)
            for k, p in DATA.items() if _HAVE[k]}


@pytest.fixture
def win(qapp):
    from pycsamt.app.desktop.windows.airborne_window import AirborneWindow

    w = AirborneWindow(parent=None)
    w.show()
    yield w
    w.close()
    plt.close("all")


# ── studio (Qt-free) ──────────────────────────────────────────────────────
class TestStudio:
    @pytest.mark.parametrize("tech", list(DATA))
    def test_every_view_runs(self, surveys, tech):
        if tech not in surveys:
            pytest.skip(f"{tech} data missing")
        s = surveys[tech]
        ctx = st.data_context(s)
        assert tech in ctx["techs"] and ctx["freqs"]
        geo = defaults(list(st.GEOMETRY_FIELDS))
        views = st.views_for(ctx["techs"])
        assert views
        for v in views:
            out = st.run_view(v, s, defaults(st.options_for(v, ctx)), geo)
            if v.kind == "table":
                assert out and all(len(df) for _t, df in out), v.key
            else:
                assert figure_blank_reason(out) is None, v.key
            plt.close("all")

    def test_views_follow_technology(self):
        ztem = {v.key for v in st.views_for(["ztem"])}
        mmt = {v.key for v in st.views_for(["mobilemt"])}
        assert ztem and mmt and ztem != mmt
        assert st.views_for([]) == [] or all(
            not v.techs for v in st.views_for([]))

    @_need("ztem")
    def test_frequency_choices_come_from_data(self, surveys):
        ctx = st.data_context(surveys["ztem"])
        for v in st.views_for(ctx["techs"]):
            f = next((f for f in st.options_for(v, ctx)
                      if f.key == "frequency_hz"), None)
            if f is not None:
                values = [c for c, _l in f.choices if c not in ("", None)]
                assert values and all(
                    any(abs(float(c) - x) < 1e-6 * max(1, x)
                        for x in ctx["freqs"]) for c in values)
                return
        pytest.skip("no frequency option for ZTEM")

    @_need("ztem")
    def test_collinear_ztem_map_draws_points(self, surveys):
        """Forrestania's stations lie on one line: the map used to crash in
        Qhull (griddata cannot triangulate collinear points)."""
        import numpy as np

        from pycsamt.emtools.ztem import plot_ztem_map

        s = surveys["ztem"]
        xy = np.array([x.coords[:2] for x in s])
        assert np.linalg.matrix_rank(xy - xy.mean(0), tol=1e-9) == 1
        ax = plot_ztem_map(s)
        fig = ax.figure if hasattr(ax, "figure") else ax
        assert figure_blank_reason(fig) is None
        plt.close("all")


# ── window ────────────────────────────────────────────────────────────────
class TestWindow:
    def test_empty_state(self, win):
        assert not win._ctrl.has_data
        assert win._view_list.count() == 0
        assert "No airborne" in win._data_status.text()

    @_need("ztem")
    def test_load_ztem(self, win):
        n = win.load(str(DATA["ztem"]))
        assert n > 0
        assert "ZTEM" in win._data_status.text()
        keys = [win._view_list.item(r).data(Qt_UserRole())
                for r in range(win._view_list.count())]
        assert set(k for k in keys if k) == {
            v.key for v in st.views_for(["ztem"])}
        assert win._stations.rowCount() == n
        assert win._status.text().startswith("✓")

    @_need("ztem")
    def test_line_filter(self, win):
        win.load(str(DATA["ztem"]))
        lines = win._ctx["lines"]
        if len(lines) < 2:
            pytest.skip("one line only")
        win._line_combo.setCurrentIndex(1)
        assert all(x.line_id == lines[0] for x in win._asites())

    @_need("ztem")
    def test_every_view_draws_in_window(self, win):
        win.load(str(DATA["ztem"]))
        for v in st.views_for(["ztem"]):
            win.select_view(v.key)
            win._on_draw()
            assert win._status.text().startswith("✓") or (
                v.kind == "table" and "table" in win._status.text()), (
                v.key, win._status.text())

    @_need("afmag_original")
    def test_geometry_group_follows_view(self, win):
        win.load(str(DATA["afmag_original"]))
        views = st.views_for(win._ctx["techs"])
        for v in views:
            win.select_view(v.key)
            assert win._grp_geo.isVisibleTo(win) == st.needs_geometry(v)

    @_need("mobilemt")
    def test_table_view_and_export(self, win, tmp_path):
        win.load(str(DATA["mobilemt"]))
        tables = [v for v in st.views_for(win._ctx["techs"])
                  if v.kind == "table"]
        if not tables:
            pytest.skip("no table view for MobileMT")
        win.select_view(tables[0].key)
        win._on_draw()
        assert win._tables and win._table.rowCount() > 0
        path = win.export_table(str(tmp_path / "t.csv"))
        assert Path(path).read_text().strip()

    @_need("ztem")
    def test_clear(self, win):
        win.load(str(DATA["ztem"]))
        win._on_clear()
        assert not win._ctrl.has_data and win._view_list.count() == 0

    def test_load_failure_reported(self, win, tmp_path):
        win._load_reporting(str(tmp_path))
        assert ("failed" in win._data_status.text().lower()
                or not win._ctrl.has_data)

    def test_dialog_cancel_is_noop(self, win, monkeypatch):
        from PySide6.QtWidgets import QFileDialog

        monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                            staticmethod(lambda *a, **k: ""))
        win._on_load()
        assert not win._ctrl.has_data


def Qt_UserRole():
    from PySide6.QtCore import Qt

    return Qt.ItemDataRole.UserRole
