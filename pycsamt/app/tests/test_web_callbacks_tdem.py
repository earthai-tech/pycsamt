# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.tdem (TDEM analysis page callbacks).

Real data
---------
data/TEMAVG/JIANGSU/ -- a small 3-station subset (TEM100, TEM1020,
TEM1060) is copied into a session tmp dir and loaded into the module's
singleton ``_CTRL`` (a real ``TDEMController``), mirroring the strategy
in ``test_tdem_controller.py``. Real ``pycsamt.tdem`` plot classes are
exercised end to end for ``run_tdem`` rather than mocked.

The one thing that must never actually run is ``tkinter.filedialog
.askdirectory`` in ``browse_tdem_folder`` -- it opens a real modal
dialog and blocks forever under a headless test run, so ``tkinter`` is
faked via ``sys.modules`` for those tests instead of exercised for real.
"""

from __future__ import annotations

import shutil
import sys
import types
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
from dash import no_update
from dash.exceptions import PreventUpdate

import pycsamt.app.web.callbacks.tdem as tdem_mod
from pycsamt.app.desktop.controllers.tdem_controller import TDEMController
from pycsamt.app.web.layout import IDs

# ── Dash callback-lookup helpers (shared pattern across web-callback tests) ──


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_multi(web_app, *substrings):
    key = next(
        k for k in web_app.callback_map if all(s in k for s in substrings)
    )
    return _unwrap(web_app.callback_map[key])


def _cb_by_input(web_app, output_substr, input_id):
    for k, v in web_app.callback_map.items():
        if output_substr not in k:
            continue
        if any(input_id in str(i.get("id")) for i in v.get("inputs", [])):
            return _unwrap(v)
    raise AssertionError(
        f"no callback found for output~={output_substr!r} input={input_id!r}"
    )


# ── Real TDEM data fixtures ───────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_JIANGSU = _ROOT / "data" / "TEMAVG" / "JIANGSU"
_HAS_JIANGSU = _JIANGSU.exists() and any(_JIANGSU.glob("*.AVG"))
_STEMS = ["TEM100", "TEM1020", "TEM1060"]


@pytest.fixture(scope="session")
def tdem_dir(tmp_path_factory):
    """A small real TDEM survey folder (3 stations) copied out of JIANGSU."""
    if not _HAS_JIANGSU:
        pytest.skip("data/TEMAVG/JIANGSU not available")
    dest = tmp_path_factory.mktemp("tdem_cb_subset")
    for stem in _STEMS:
        for f in _JIANGSU.glob(stem + ".*"):
            shutil.copy(f, dest / f.name)
    return dest


@pytest.fixture(scope="session")
def loaded_ctrl_state(tdem_dir):
    """
    Load real data into the module-level ``_CTRL`` singleton once, and
    hand back the folder path for tests that need it. Loading is real
    end-to-end (``read_temavg_survey`` / ``read_temavg_soundings``).
    """
    ok = tdem_mod._CTRL.load_folder(str(tdem_dir))
    assert ok is True
    yield str(tdem_dir)


def _fig():
    return plt.figure(figsize=(6, 4))


def _close():
    plt.close("all")


class _FakeCtx:
    def __init__(self, triggered_id):
        self.triggered_id = triggered_id


# ── 1. browse_tdem_folder ─────────────────────────────────────────────────────


class _FakeTkRoot:
    def __init__(self, scaling_raises=False, attrs_raise=False):
        self._scaling_raises = scaling_raises
        self._attrs_raise = attrs_raise
        self.tk = types.SimpleNamespace(call=self._tk_call)
        self.destroyed = False

    def _tk_call(self, *_a, **_k):
        if self._scaling_raises:
            raise RuntimeError("no scaling support")

    def withdraw(self):
        pass

    def attributes(self, *_a, **_k):
        if self._attrs_raise:
            raise RuntimeError("no topmost support")

    def update(self):
        pass

    def destroy(self):
        self.destroyed = True


def _install_fake_tkinter(monkeypatch, *, folder="", root=None, tk_raises=False):
    fake_tk = types.ModuleType("tkinter")
    fake_filedialog = types.ModuleType("tkinter.filedialog")
    fake_font = types.ModuleType("tkinter.font")

    if tk_raises:

        def _tk_ctor():
            raise RuntimeError("no display available")

    else:
        _root = root if root is not None else _FakeTkRoot()

        def _tk_ctor():
            return _root

    fake_tk.Tk = _tk_ctor
    fake_filedialog.askdirectory = lambda **_kw: folder
    fake_font.nametofont = lambda name: types.SimpleNamespace(
        cget=lambda k: 9, configure=lambda **kw: None
    )
    fake_tk.filedialog = fake_filedialog
    fake_tk.font = fake_font

    monkeypatch.setitem(sys.modules, "tkinter", fake_tk)
    monkeypatch.setitem(sys.modules, "tkinter.filedialog", fake_filedialog)
    monkeypatch.setitem(sys.modules, "tkinter.font", fake_font)


class TestBrowseTdemFolder:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.TDEM_FOLDER}.value")

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None)

    def test_dialog_returns_folder(self, monkeypatch, web_app):
        _install_fake_tkinter(monkeypatch, folder="D:/some/tdem/folder")
        out = self._fn(web_app)(1)
        assert out == "D:/some/tdem/folder"

    def test_dialog_cancelled_returns_no_update(self, monkeypatch, web_app):
        _install_fake_tkinter(monkeypatch, folder="")
        out = self._fn(web_app)(1)
        assert out is no_update

    def test_scaling_failure_is_swallowed(self, monkeypatch, web_app):
        # root.tk.call("tk", "scaling", ...) has its own try/except pass,
        # so a failure there must not stop the dialog flow.
        root = _FakeTkRoot(scaling_raises=True)
        _install_fake_tkinter(monkeypatch, folder="X:/folder", root=root)
        out = self._fn(web_app)(1)
        assert out == "X:/folder"
        assert root.destroyed is True

    def test_attributes_failure_is_not_swallowed(self, monkeypatch, web_app):
        # Unlike the scaling call, root.attributes("-topmost", ...) is NOT
        # wrapped in its own try/except -- a failure there propagates to
        # the function's outer except, which converts it to PreventUpdate
        # (root.destroy() is then never reached).
        root = _FakeTkRoot(attrs_raise=True)
        _install_fake_tkinter(monkeypatch, folder="X:/folder", root=root)
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1)
        assert root.destroyed is False

    def test_tk_unavailable_prevents_update(self, monkeypatch, web_app):
        _install_fake_tkinter(monkeypatch, tk_raises=True)
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1)


# ── 2. switch_tdem_tab ────────────────────────────────────────────────────────


class TestSwitchTdemTab:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.TDEM_ACTIVE_TAB}.data")

    def test_no_trigger_no_update(self, monkeypatch, web_app):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(None))
        assert self._fn(web_app)(1, None, None, None) is no_update

    def test_decay_tab_clicked(self, monkeypatch, web_app):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx("tdem-tab-decay"))
        assert self._fn(web_app)(1, None, None, None) == "decay"

    def test_dashboard_tab_clicked(self, monkeypatch, web_app):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx("tdem-tab-dashboard"))
        assert self._fn(web_app)(None, None, None, 1) == "dashboard"


# ── 3. sync_tdem_plots ────────────────────────────────────────────────────────


class TestSyncTdemPlots:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.TDEM_PLOT}.options", f"{IDs.TDEM_PLOT}.value")

    def test_no_active_tab_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None)

    def test_decay_tab_options(self, web_app):
        opts, value = self._fn(web_app)("decay")
        assert {"label": "Decay curves", "value": "PlotDecayCurve"} in opts
        assert value == opts[0]["value"]

    def test_section_tab_options(self, web_app):
        opts, value = self._fn(web_app)("section")
        values = {o["value"] for o in opts}
        assert "PlotTEMAVGSection" in values
        assert value == opts[0]["value"]

    def test_unknown_tab_falls_back_to_decay_plots(self, web_app):
        opts, _value = self._fn(web_app)("not-a-real-tab")
        values = {o["value"] for o in opts}
        assert values == {"PlotDecayCurve", "PlotTransformedRho"}


# ── 4. update_param_panels ────────────────────────────────────────────────────


class TestUpdateParamPanels:
    def _fn(self, web_app):
        return _cb_multi(
            web_app,
            "tdem-sounding-section.style",
            "tdem-cmap-section.style",
            "tdem-gate-section.style",
        )

    def test_no_class_name_all_hidden(self, web_app):
        sounding, cmap, gate = self._fn(web_app)(None)
        assert sounding == cmap == gate == {"display": "none"}

    def test_sounding_plot_shows_sounding_panel_only(self, web_app):
        sounding, cmap, gate = self._fn(web_app)("PlotDecayCurve")
        assert sounding == {"display": "block"}
        assert cmap == {"display": "none"}
        assert gate == {"display": "none"}

    def test_section_cmap_plot_shows_cmap_panel(self, web_app):
        sounding, cmap, gate = self._fn(web_app)("PlotTEMAVGSection")
        assert sounding == {"display": "none"}
        assert cmap == {"display": "block"}
        assert gate == {"display": "none"}

    def test_gate_profile_shows_both_cmap_and_gate_panels(self, web_app):
        # PlotGateProfile is in both _SECTION_CMAP_PLOTS and _GATE_PLOTS.
        sounding, cmap, gate = self._fn(web_app)("PlotGateProfile")
        assert sounding == {"display": "none"}
        assert cmap == {"display": "block"}
        assert gate == {"display": "block"}


# ── 5. update_ctx_bar ─────────────────────────────────────────────────────────


class TestUpdateCtxBar:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.TDEM_CTX_BAR}.children")

    def test_defaults_to_decay_group_no_plot_no_data(self, web_app):
        parts = self._fn(web_app)(None, None, None)
        group_span = parts[0]
        assert "Decay / Rho" in group_span.children
        # No plot label -> only the group span + the "not loaded" hint.
        assert len(parts) == 2
        assert "Browse and load" in parts[1].children

    def test_with_plot_label_adds_separator_and_label(self, web_app):
        parts = self._fn(web_app)("decay", "PlotDecayCurve", None)
        texts = [
            p.children if isinstance(p.children, str) else None for p in parts
        ]
        assert "Decay curves" in texts
        assert "Browse and load" in parts[-1].children

    def test_data_loaded_shows_loaded_stat(self, web_app):
        parts = self._fn(web_app)("map", "PlotSurveyMap", True)
        joined = " ".join(
            p.children for p in parts if isinstance(p.children, str)
        )
        assert "data loaded" in joined

    def test_unknown_tab_falls_back_to_decay_group(self, web_app):
        parts = self._fn(web_app)("not-a-tab", None, None)
        assert "Decay / Rho" in parts[0].children


# ── 6. load_tdem ──────────────────────────────────────────────────────────────


class TestLoadTdem:
    def _fn(self, web_app):
        return _cb_multi(
            web_app, f"{IDs.TDEM_INFO}.children", f"{IDs.STORE_TDEM_LOADED}.data"
        )

    def test_no_folder_no_sample_prevents_update(self, monkeypatch, web_app):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_LOAD))
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, None, None)

    def test_sample_button_uses_sample_path(self, monkeypatch, web_app, tdem_dir):
        # Point the module's sample-data constant at our small real subset
        # instead of the full 55-station JIANGSU folder so this stays fast.
        monkeypatch.setattr(tdem_mod, "_SAMPLE", str(tdem_dir))
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_SAMPLE))
        info, loaded, folder, is_open, body = self._fn(web_app)(None, 1, None)
        assert loaded is True
        assert folder == str(tdem_dir)
        assert is_open is False
        assert body == ""
        assert len(info) > 0  # _fmt_info() produced info-line Divs

    def test_real_folder_load_success(self, monkeypatch, web_app, tdem_dir):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_LOAD))
        info, loaded, folder, is_open, body = self._fn(web_app)(
            1, None, str(tdem_dir) + "  "
        )
        assert loaded is True
        assert folder == str(tdem_dir)  # stripped
        assert is_open is False

    def test_empty_folder_reports_no_files_found(self, monkeypatch, web_app, tmp_path):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_LOAD))
        info, loaded, folder, is_open, body = self._fn(web_app)(
            1, None, str(tmp_path)
        )
        assert loaded is False
        assert info == "No TDEM files found in that folder."
        assert is_open is True
        assert "TDEM: folder empty" in body

    def test_load_folder_exception_is_caught(self, monkeypatch, web_app):
        # TDEMController.load_folder() never raises in practice (it wraps
        # its own internals in try/except and always returns bool), so the
        # outer try/except in load_tdem() is otherwise unreachable through
        # real data -- force it directly to exercise that defensive branch.
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_LOAD))

        def _boom(_folder):
            raise RuntimeError("disk exploded")

        monkeypatch.setattr(tdem_mod._CTRL, "load_folder", _boom)
        info, loaded, folder, is_open, body = self._fn(web_app)(
            1, None, "Z:/wherever"
        )
        assert loaded is False
        assert "Error: disk exploded" in info
        assert is_open is True
        assert "TDEM load: disk exploded" in body


# ── 7. run_tdem ────────────────────────────────────────────────────────────────


class TestRunTdem:
    def _fn(self, web_app):
        return _cb_multi(
            web_app, f"{IDs.IMG_TDEM_DECAY}.src", f"{IDs.TDEM_SPINNER}.children"
        )

    def _base_kwargs(self):
        return dict(
            n_clicks=1,
            class_name="PlotDecayCurve",
            nav_section="tdem",
            active_tab="decay",
            tdem_loaded=True,
            theme="dark",
            figsize="10x5",
            cmap="viridis",
            sounding_idx=0,
            gate_idx=0,
        )

    def test_plot_trigger_off_tdem_section_no_update(self, monkeypatch, web_app):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.TDEM_PLOT))
        kw = self._base_kwargs()
        kw["nav_section"] = "home"
        out = self._fn(web_app)(**kw)
        assert out == (no_update,) * 7

    def test_nav_section_trigger_not_loaded_no_update(self, monkeypatch, web_app):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.NAV_SECTION))
        kw = self._base_kwargs()
        kw["tdem_loaded"] = False
        out = self._fn(web_app)(**kw)
        assert out == (no_update,) * 7

    def test_no_class_name_prevents_update(self, monkeypatch, web_app):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_RUN))
        kw = self._base_kwargs()
        kw["class_name"] = None
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(**kw)

    def test_not_loaded_ctrl_shows_error_toast(self, monkeypatch, web_app):
        # Fresh, never-loaded controller instance -- swap it in temporarily.
        monkeypatch.setattr(tdem_mod, "_CTRL", TDEMController())
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_RUN))
        kw = self._base_kwargs()
        out = self._fn(web_app)(**kw)
        decay_src, section_src, map_src, dash_src, spinner, is_open, body = out
        assert decay_src == section_src == map_src == dash_src
        assert spinner == ""
        assert is_open is True
        assert "load a folder first" in body

    def test_successful_decay_render_dark(
        self, monkeypatch, web_app, loaded_ctrl_state
    ):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_RUN))
        kw = self._base_kwargs()
        out = self._fn(web_app)(**kw)
        decay_src, section_src, map_src, dash_src, spinner, is_open, body = out
        print("DEBUG body=", repr(body), "is_open=", is_open, "decay_src=", repr(decay_src)[:80])
        assert decay_src.startswith("data:image/png;base64,")
        assert section_src is no_update
        assert map_src is no_update
        assert dash_src is no_update
        assert spinner == ""
        assert is_open is False
        assert body == ""

    def test_successful_render_light_theme_routes_to_active_tab_slot(
        self, monkeypatch, web_app, loaded_ctrl_state
    ):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_RUN))
        kw = self._base_kwargs()
        kw["theme"] = "light"
        kw["class_name"] = "PlotTEMAVGSection"
        kw["active_tab"] = "section"
        out = self._fn(web_app)(**kw)
        decay_src, section_src, map_src, dash_src, *_ = out
        assert decay_src is no_update
        assert section_src.startswith("data:image/png;base64,")
        assert map_src is no_update
        assert dash_src is no_update

    def test_gate_profile_dashboard_index_slot(
        self, monkeypatch, web_app, loaded_ctrl_state
    ):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_RUN))
        kw = self._base_kwargs()
        kw["class_name"] = "PlotGateProfile"
        kw["active_tab"] = "section"
        kw["gate_idx"] = 1
        out = self._fn(web_app)(**kw)
        assert out[1].startswith("data:image/png;base64,")

    def test_dashboard_render(self, monkeypatch, web_app, loaded_ctrl_state):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_RUN))
        kw = self._base_kwargs()
        kw["class_name"] = "PlotTEMDashboard"
        kw["active_tab"] = "dashboard"
        out = self._fn(web_app)(**kw)
        assert out[3].startswith("data:image/png;base64,")

    def test_bad_figsize_falls_back_to_default(
        self, monkeypatch, web_app, loaded_ctrl_state
    ):
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_RUN))
        kw = self._base_kwargs()
        kw["figsize"] = "not-a-real-size"
        out = self._fn(web_app)(**kw)
        assert out[0].startswith("data:image/png;base64,")

    def test_invalid_sounding_idx_hits_exception_branch(
        self, monkeypatch, web_app, loaded_ctrl_state
    ):
        # int("abc") raises ValueError while building the sounding_index
        # kwarg, which is caught by run_tdem()'s outer except -> error toast.
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_RUN))
        kw = self._base_kwargs()
        kw["sounding_idx"] = "abc"
        out = self._fn(web_app)(**kw)
        decay_src, section_src, map_src, dash_src, spinner, is_open, body = out
        assert decay_src is no_update
        assert spinner == ""
        assert is_open is True
        assert "PlotDecayCurve:" in body

    def test_unloaded_class_name_still_renders_via_placeholder_annotation(
        self, monkeypatch, web_app, loaded_ctrl_state
    ):
        # TDEMController.draw() never raises for an unknown class name --
        # it annotates the figure with an error message and returns None,
        # so run_tdem() still produces a normal (non-error) image.
        monkeypatch.setattr(tdem_mod, "ctx", _FakeCtx(IDs.BTN_TDEM_RUN))
        kw = self._base_kwargs()
        kw["class_name"] = "NotARealPlotClass"
        out = self._fn(web_app)(**kw)
        assert out[0].startswith("data:image/png;base64,")
        assert out[5] is False


# ── _ctrl_draw ────────────────────────────────────────────────────────────────


class TestCtrlDraw:
    def test_no_kwargs_calls_draw_directly(self, monkeypatch):
        calls = []
        ctrl = TDEMController()
        monkeypatch.setattr(
            tdem_mod,
            "_CTRL",
            types.SimpleNamespace(
                draw=lambda *a, **k: calls.append(("draw", a, k)) or "ok"
            ),
        )
        fig = _fig()
        result = tdem_mod._ctrl_draw("PlotDecayCurve", True, "soundings", fig)
        assert result == "ok"
        assert calls == [("draw", ("PlotDecayCurve", True, "soundings", fig), {})]
        _close()

    def test_kwargs_without_draw_ext_silently_dropped(self, monkeypatch):
        """
        BUG (documented, not fixed): the real TDEMController has no
        ``draw_ext`` method, so ``hasattr(_CTRL, "draw_ext")`` is always
        False in production. Every extra parameter the TDEM page collects
        -- cmap, sounding_index, gate_index -- is therefore silently
        dropped on every real render; only the plain ``draw()`` 4-arg
        signature is ever actually invoked, regardless of what kwargs the
        caller passes in.
        """
        assert not hasattr(tdem_mod._CTRL, "draw_ext")
        calls = []
        fake_ctrl = types.SimpleNamespace(
            draw=lambda *a, **k: calls.append(("draw", a, k)) or "ok"
        )
        monkeypatch.setattr(tdem_mod, "_CTRL", fake_ctrl)
        fig = _fig()
        result = tdem_mod._ctrl_draw(
            "PlotTEMAVGSection", True, "survey", fig, cmap="plasma"
        )
        assert result == "ok"
        # draw_ext was never called (fake controller doesn't even have one);
        # draw() was called instead, and the cmap kwarg never reached it.
        assert calls == [
            ("draw", ("PlotTEMAVGSection", True, "survey", fig), {})
        ]
        _close()

    def test_draw_ext_used_when_present(self, monkeypatch):
        calls = []

        class _Ctrl:
            def draw_ext(self, *a, **k):
                calls.append(("draw_ext", a, k))
                return "ext-ok"

            def draw(self, *a, **k):
                calls.append(("draw", a, k))
                return "plain-ok"

        monkeypatch.setattr(tdem_mod, "_CTRL", _Ctrl())
        fig = _fig()
        result = tdem_mod._ctrl_draw(
            "PlotDecayCurve", True, "soundings", fig, sounding_index=2
        )
        assert result == "ext-ok"
        assert calls == [
            (
                "draw_ext",
                ("PlotDecayCurve", True, "soundings", fig),
                {"sounding_index": 2},
            )
        ]
        _close()

    def test_draw_ext_type_error_falls_back_to_draw(self, monkeypatch):
        calls = []

        class _Ctrl:
            def draw_ext(self, *a, **k):
                calls.append("draw_ext")
                raise TypeError("unexpected kwarg")

            def draw(self, *a, **k):
                calls.append("draw")
                return "fallback-ok"

        monkeypatch.setattr(tdem_mod, "_CTRL", _Ctrl())
        fig = _fig()
        result = tdem_mod._ctrl_draw(
            "PlotDecayCurve", True, "soundings", fig, sounding_index=2
        )
        assert result == "fallback-ok"
        assert calls == ["draw_ext", "draw"]
        _close()


# ── _fmt_info ─────────────────────────────────────────────────────────────────


class TestFmtInfo:
    def test_empty_summary_returns_empty_list(self):
        assert tdem_mod._fmt_info("") == []

    def test_none_summary_returns_empty_list(self):
        assert tdem_mod._fmt_info(None) == []

    def test_multiline_summary_produces_one_div_per_nonempty_line(self):
        summary = "Folder: X\nAVG files: 3   Z files: 3\nSoundings: 3\n"
        divs = tdem_mod._fmt_info(summary)
        assert len(divs) == 3
        assert divs[0].children == "Folder: X"
        assert all(d.className == "tdem-info-line" for d in divs)

    def test_blank_lines_are_filtered_out(self):
        summary = "Line one\n\n\nLine two\n"
        divs = tdem_mod._fmt_info(summary)
        assert [d.children for d in divs] == ["Line one", "Line two"]
