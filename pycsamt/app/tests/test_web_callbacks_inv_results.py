# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for pycsamt.app.web.callbacks.inv_results (Inversion Results Viewer).

Strategy
--------
Real, fast, shipped sample inversion outputs are used instead of mocks:
  * Occam2D  -> ``data/occam2D/``            (real InversionResult, ~1s load)
  * MARE2DEM -> ``data/mare2dem/demo_mt_inversion/`` (~0.1s load)
ModEM needs an external ``ModEMv626`` example tree that isn't shipped with
this repo (same skip convention as ``pycsamt/models/modem/tests``), so
ModEM-specific plotting/loading is skipped here when that data is absent;
the solver-agnostic parts of the callback layer (folder browser, generic
dispatch/caching) are still exercised without it.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
from dash.exceptions import PreventUpdate

import pycsamt.app.web.callbacks.inv_results as ivr
from pycsamt.app.web.layout import IDs

_ROOT = Path(__file__).parents[3]
_OCCAM_DIR = _ROOT / "data" / "occam2D"
_MARE_DIR = _ROOT / "data" / "mare2dem" / "demo_mt_inversion"
_MODEM_DIR = _ROOT / "ModEMv626" / "ModEM" / "examples" / "2D_MT" / "BLOCK2"

_HAS_OCCAM = _OCCAM_DIR.exists()
_HAS_MARE = _MARE_DIR.exists()
_HAS_MODEM = _MODEM_DIR.exists()


def _default_ctrl_args(**overrides):
    """32 positional values matching _control_states()/_build_controls order."""
    defaults = dict(
        which="final",
        iteration=None,
        rho_min=1.0,
        rho_max=1000.0,
        depth_max=2000.0,
        cmap="jet_r",
        sect_mode="axial",
        sect_dir="NS",
        sect_offset=0,
        sect_start_lat=None,
        sect_start_lon=None,
        sect_end_lat=None,
        sect_end_lon=None,
        sect_n_samp=200,
        sect_orig_lat=None,
        sect_orig_lon=None,
        dm_depths="50,200,500,1000",
        dm_ncols=2,
        dm_lat=None,
        dm_lon=None,
        ap_dir="NS",
        ap_n=5,
        ap_offsets="",
        ap_names="",
        ap_ncols=3,
        cov_cmap="Blues",
        resp_station=None,
        resp_comp="ZXY",
        ps_comp="ZXY",
        show_sta=True,
        show_names=True,
        sta_tol=300.0,
    )
    defaults.update(overrides)
    return tuple(defaults.values())


@pytest.fixture(scope="session")
def occam_result():
    if not _HAS_OCCAM:
        pytest.skip("occam2D sample data not available")
    from pycsamt.models.occam2d.results import InversionResult

    return InversionResult(workdir=_OCCAM_DIR)


@pytest.fixture(scope="session")
def mare_result():
    if not _HAS_MARE:
        pytest.skip("MARE2DEM sample data not available")
    from pycsamt.models.mare2dem.results import InversionResult

    return InversionResult(workdir=_MARE_DIR)


@pytest.fixture(autouse=True)
def _clear_result_cache():
    yield
    ivr._RESULT_CACHE.clear()


# ── Fake ModEM result + plot classes ────────────────────────────────────────
# ModEM needs an external ModEMv626 example tree that isn't shipped with this
# repo (same skip convention as pycsamt/models/modem/tests), so the
# ModEM-specific *branch logic* inside _make_fig_modem/_load_result -- model
# selection, control parsing, arg building, dispatch -- is instead exercised
# against fake result objects and stubbed pycsamt.models.modem.plot classes.
# This tests inv_results.py's own logic, not the real plotting code.


class _FakeModel:
    def __init__(self, nx=8, ny=8, nz=4):
        self.nx, self.ny, self.nz = nx, ny, nz


class _FakeModemResult:
    def __init__(self):
        self.model_final = _FakeModel()
        self.model_initial = _FakeModel(4, 4, 4)
        self.models = {
            "iter_001": _FakeModel(1, 1, 1),
            "iter_002": _FakeModel(2, 2, 2),
            "iter_010": _FakeModel(10, 10, 10),
        }
        self.mode = "TE_TM"
        self.n_iter = 10
        self.final_rms = 1.532
        self.data_obs = SimpleNamespace(stations=["mt02", "mt01"])
        self.data_pred = None


class _FakePlot:
    """Records constructor kwargs; ``.plot()`` returns a throwaway Figure."""

    last_kwargs = None

    def __init__(self, *args, **kwargs):
        type(self).last_kwargs = kwargs
        self.args = args
        self.kwargs = kwargs

    def plot(self, **kw):
        import matplotlib.pyplot as plt

        fig, _ = plt.subplots()
        return fig


@pytest.fixture
def fake_modem_plots(monkeypatch):
    import pycsamt.models.modem.plot as modem_plot

    names = [
        "PlotAllProfiles",
        "PlotCovariance",
        "PlotDepthMap",
        "PlotMisfit",
        "PlotPseudo",
        "PlotResponse",
        "PlotSection",
    ]
    fakes = {}
    for name in names:
        fake_cls = type(name, (_FakePlot,), {"last_kwargs": None})
        fakes[name] = fake_cls
        monkeypatch.setattr(modem_plot, name, fake_cls)
    return fakes


@pytest.fixture
def fake_modem_result():
    return _FakeModemResult()


def _unwrap(entry):
    fn = entry["callback"]
    return getattr(fn, "__wrapped__", fn)


def _cb(web_app, output_id_prop):
    return _unwrap(web_app.callback_map[output_id_prop])


def _cb_multi(web_app, *substrings):
    key = next(k for k in web_app.callback_map if all(s in k for s in substrings))
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


# ── Pure helpers ─────────────────────────────────────────────────────────────


class TestListDirs:
    def test_lists_subfolders_sorted(self, tmp_path):
        (tmp_path / "b").mkdir()
        (tmp_path / "a").mkdir()
        (tmp_path / "file.txt").write_text("x")
        (tmp_path / ".hidden").mkdir()
        assert ivr._list_dirs(str(tmp_path)) == ["a", "b"]

    def test_nonexistent_path_returns_empty(self):
        assert ivr._list_dirs("/no/such/path/xyz") == []


class TestDirListing:
    def test_empty_folder_shows_hint(self, tmp_path):
        items = ivr._dir_listing(str(tmp_path))
        assert "No sub-folders" in items[0].children

    def test_lists_items(self, tmp_path):
        (tmp_path / "sub1").mkdir()
        items = ivr._dir_listing(str(tmp_path))
        assert len(items) == 1  # one ListGroup wrapping all rows


class TestDetectSolver:
    def test_not_a_dir(self, tmp_path):
        assert ivr._detect_solver(tmp_path / "nope") == "unknown"

    def test_modem_by_rho_file(self, tmp_path):
        (tmp_path / "model.rho").write_text("x")
        assert ivr._detect_solver(tmp_path) == "modem"

    def test_occam2d_by_startup_file(self, tmp_path):
        (tmp_path / "startup").write_text("x")
        assert ivr._detect_solver(tmp_path) == "occam2d"

    def test_mare2dem_by_emdata_file(self, tmp_path):
        (tmp_path / "run.emdata").write_text("x")
        assert ivr._detect_solver(tmp_path) == "mare2dem"

    def test_unknown_empty_dir(self, tmp_path):
        assert ivr._detect_solver(tmp_path) == "unknown"

    @pytest.mark.skipif(not _HAS_OCCAM, reason="occam2D data not available")
    def test_real_occam2d_dir(self):
        assert ivr._detect_solver(_OCCAM_DIR) == "occam2d"

    @pytest.mark.skipif(not _HAS_MARE, reason="MARE2DEM data not available")
    def test_real_mare2dem_dir(self):
        assert ivr._detect_solver(_MARE_DIR) == "mare2dem"


class TestResultCache:
    def test_store_and_fetch(self):
        ivr._store_result("sess-1", {"x": 1})
        assert ivr._fetch_result("sess-1") == {"x": 1}

    def test_fetch_missing_returns_none(self):
        assert ivr._fetch_result("no-such-session") is None


class TestBuildControls:
    def test_maps_positional_args_to_dict(self):
        args = _default_ctrl_args(resp_station="18-001A", show_sta=False)
        controls = ivr._build_controls(*args)
        assert controls["resp_station"] == "18-001A"
        assert controls["show_sta"] is False
        assert controls["cmap"] == "jet_r"


class TestBuildInfoStrip:
    def test_occam2d(self, occam_result):
        items = ivr._build_info_strip(occam_result, "occam2d", str(_OCCAM_DIR))
        texts = [str(getattr(i, "children", i)) for i in items]
        assert any("OCCAM2D" in t for t in texts)
        assert any("RMS=" in t for t in texts)

    def test_mare2dem(self, mare_result):
        items = ivr._build_info_strip(mare_result, "mare2dem", str(_MARE_DIR))
        texts = [str(getattr(i, "children", i)) for i in items]
        assert any("MARE2DEM" in t for t in texts)

    def test_unknown_solver_swallows_and_returns_base_items(self):
        items = ivr._build_info_strip(object(), "totally-unknown", "/tmp/x")
        assert len(items) == 2  # solver badge + path, no RMS/iter badges

    def test_modem(self):
        result = _FakeModemResult()
        items = ivr._build_info_strip(result, "modem", "/tmp/BLOCK2")
        texts = [str(getattr(i, "children", i)) for i in items]
        assert any("MODEM" in t for t in texts)
        assert any("grid=8" in t and "8" in t and "4" in t for t in texts)
        assert any("iter=10" in t for t in texts)
        assert any("RMS=1.532" in t for t in texts)

    def test_modem_badge_exception_is_swallowed(self):
        result = _FakeModemResult()
        result.final_rms = "not-a-number"  # f"{rms:.3f}" raises TypeError
        items = ivr._build_info_strip(result, "modem", "/tmp/x")
        assert len(items) == 2


# ── Real plot dispatch (Occam2D + MARE2DEM) ─────────────────────────────────


@pytest.mark.skipif(not _HAS_OCCAM, reason="occam2D data not available")
class TestMakeFigOccam2d:
    def test_all_tabs_render(self, occam_result):
        controls = ivr._build_controls(*_default_ctrl_args())
        for tab in [
            "conv",
            "model",
            "response",
            "pseudo",
            "sounding1d",
            "sitemisfit",
            "respgrid",
        ]:
            fig = ivr._make_fig_occam2d(occam_result, tab, controls)
            assert fig is not None

    def test_unknown_tab_raises(self, occam_result):
        controls = ivr._build_controls(*_default_ctrl_args())
        with pytest.raises(ValueError, match="Unknown Occam2D tab"):
            ivr._make_fig_occam2d(occam_result, "not-a-tab", controls)


class TestMakeFigModem:
    """Exercises _make_fig_modem's branch logic against fake plot classes."""

    def test_conv_tab(self, fake_modem_result, fake_modem_plots):
        controls = ivr._build_controls(*_default_ctrl_args())
        fig = ivr._make_fig_modem(fake_modem_result, "conv", controls)
        assert fig is not None
        assert fake_modem_plots["PlotMisfit"].last_kwargs["result"] is (
            fake_modem_result
        )

    def test_which_initial_selects_initial_model(
        self, fake_modem_result, fake_modem_plots
    ):
        controls = ivr._build_controls(*_default_ctrl_args(which="initial"))
        ivr._make_fig_modem(fake_modem_result, "section", controls)
        kw = fake_modem_plots["PlotSection"].last_kwargs
        assert kw["model"] is fake_modem_result.model_initial

    def test_which_iteration_selects_indexed_model(
        self, fake_modem_result, fake_modem_plots
    ):
        controls = ivr._build_controls(
            *_default_ctrl_args(which="iteration", iteration=2)
        )
        ivr._make_fig_modem(fake_modem_result, "section", controls)
        kw = fake_modem_plots["PlotSection"].last_kwargs
        assert kw["model"] is fake_modem_result.models["iter_002"]

    def test_which_iteration_out_of_range_is_clamped(
        self, fake_modem_result, fake_modem_plots
    ):
        controls = ivr._build_controls(
            *_default_ctrl_args(which="iteration", iteration=99)
        )
        ivr._make_fig_modem(fake_modem_result, "section", controls)
        kw = fake_modem_plots["PlotSection"].last_kwargs
        assert kw["model"] is fake_modem_result.models["iter_010"]

    def test_which_iteration_no_iter_keys_falls_back_to_final(
        self, fake_modem_plots
    ):
        result = _FakeModemResult()
        result.models = {}
        controls = ivr._build_controls(
            *_default_ctrl_args(which="iteration", iteration=1)
        )
        ivr._make_fig_modem(result, "section", controls)
        kw = fake_modem_plots["PlotSection"].last_kwargs
        assert kw["model"] is result.model_final

    def test_which_iteration_with_none_iteration_falls_back_to_final(
        self, fake_modem_result, fake_modem_plots
    ):
        controls = ivr._build_controls(
            *_default_ctrl_args(which="iteration", iteration=None)
        )
        ivr._make_fig_modem(fake_modem_result, "section", controls)
        kw = fake_modem_plots["PlotSection"].last_kwargs
        assert kw["model"] is fake_modem_result.model_final

    def test_model_none_falls_back_to_initial(self, fake_modem_plots):
        result = _FakeModemResult()
        result.model_final = None
        controls = ivr._build_controls(*_default_ctrl_args(which="final"))
        ivr._make_fig_modem(result, "section", controls)
        kw = fake_modem_plots["PlotSection"].last_kwargs
        assert kw["model"] is result.model_initial

    def test_section_axial(self, fake_modem_result, fake_modem_plots):
        controls = ivr._build_controls(
            *_default_ctrl_args(sect_mode="axial", sect_dir="EW", sect_offset=50)
        )
        ivr._make_fig_modem(fake_modem_result, "section", controls)
        kw = fake_modem_plots["PlotSection"].last_kwargs
        assert kw["direction"] == "EW"
        assert kw["offset"] == 50.0

    def test_section_arbitrary_complete(self, fake_modem_result, fake_modem_plots):
        controls = ivr._build_controls(
            *_default_ctrl_args(
                sect_mode="arbitrary",
                sect_start_lat=1.0,
                sect_start_lon=2.0,
                sect_end_lat=3.0,
                sect_end_lon=4.0,
                sect_orig_lat=0.5,
                sect_orig_lon=0.6,
                sect_n_samp=100,
            )
        )
        ivr._make_fig_modem(fake_modem_result, "section", controls)
        kw = fake_modem_plots["PlotSection"].last_kwargs
        assert kw["start_point"] == (1.0, 2.0)
        assert kw["end_point"] == (3.0, 4.0)
        assert kw["use_latlon"] is True
        assert kw["origin_lat"] == 0.5
        assert kw["n_samples"] == 100

    def test_section_arbitrary_missing_coords_raises(
        self, fake_modem_result, fake_modem_plots
    ):
        controls = ivr._build_controls(*_default_ctrl_args(sect_mode="arbitrary"))
        with pytest.raises(ValueError, match="Provide all four"):
            ivr._make_fig_modem(fake_modem_result, "section", controls)

    def test_depthmap(self, fake_modem_result, fake_modem_plots):
        controls = ivr._build_controls(
            *_default_ctrl_args(
                dm_depths="10,20,30", dm_ncols=3, dm_lat=1.1, dm_lon=2.2
            )
        )
        ivr._make_fig_modem(fake_modem_result, "depthmap", controls)
        kw = fake_modem_plots["PlotDepthMap"].last_kwargs
        assert kw["depths"] == [10.0, 20.0, 30.0]
        assert kw["ncols"] == 3
        assert kw["origin_lat"] == 1.1
        assert kw["origin_lon"] == 2.2

    def test_profiles_with_offsets_and_names(
        self, fake_modem_result, fake_modem_plots
    ):
        controls = ivr._build_controls(
            *_default_ctrl_args(
                ap_dir="EW", ap_n=7, ap_offsets="10, 20", ap_names="A, B", ap_ncols=4
            )
        )
        ivr._make_fig_modem(fake_modem_result, "profiles", controls)
        kw = fake_modem_plots["PlotAllProfiles"].last_kwargs
        assert kw["offsets"] == [10.0, 20.0]
        assert kw["profile_names"] == ["A", "B"]
        assert kw["n_profiles"] == 7
        assert kw["ncols"] == 4

    def test_profiles_without_offsets_and_names(
        self, fake_modem_result, fake_modem_plots
    ):
        controls = ivr._build_controls(
            *_default_ctrl_args(ap_offsets="", ap_names="")
        )
        ivr._make_fig_modem(fake_modem_result, "profiles", controls)
        kw = fake_modem_plots["PlotAllProfiles"].last_kwargs
        assert kw["offsets"] is None
        assert kw["profile_names"] is None

    def test_covariance(self, fake_modem_result, fake_modem_plots):
        controls = ivr._build_controls(*_default_ctrl_args(cov_cmap="Reds"))
        ivr._make_fig_modem(fake_modem_result, "covariance", controls)
        assert (
            fake_modem_plots["PlotCovariance"].last_kwargs["mask_cmap"] == "Reds"
        )

    def test_response(self, fake_modem_result, fake_modem_plots):
        controls = ivr._build_controls(
            *_default_ctrl_args(resp_station="mt01", resp_comp="ZXX")
        )
        ivr._make_fig_modem(fake_modem_result, "response", controls)
        kw = fake_modem_plots["PlotResponse"].last_kwargs
        assert kw["station"] == "mt01"
        assert kw["component"] == "ZXX"

    def test_pseudo(self, fake_modem_result, fake_modem_plots):
        controls = ivr._build_controls(*_default_ctrl_args(ps_comp="ZYX"))
        ivr._make_fig_modem(fake_modem_result, "pseudo", controls)
        assert fake_modem_plots["PlotPseudo"].last_kwargs["component"] == "ZYX"

    def test_unknown_tab_raises(self, fake_modem_result, fake_modem_plots):
        controls = ivr._build_controls(*_default_ctrl_args())
        with pytest.raises(ValueError, match="Unknown ModEM tab"):
            ivr._make_fig_modem(fake_modem_result, "not-a-tab", controls)


@pytest.mark.skipif(not _HAS_MARE, reason="MARE2DEM data not available")
class TestMakeFigMare2dem:
    def test_conv_model_response_render(self, mare_result):
        for tab in ["conv", "model", "response"]:
            fig = ivr._make_fig_mare2dem(mare_result, tab, {})
            assert fig is not None

    def test_survey_tab_is_broken_for_plain_inversion_result(self, mare_result):
        """Real bug: the "survey" tab does
        ``em_arg = result.em if hasattr(result, "em") else result`` then
        ``PlotSurveyLayout(em_arg).plot(...)``. A plain
        ``mare2dem.results.InversionResult`` has neither an ``.em``
        attribute nor the ``.utm`` attribute ``PlotSurveyLayout`` needs, so
        this always raises -- the web callback's own try/except swallows
        it into a generic "Plot error" placeholder figure (see
        ``_generate_plot``), but the survey-layout tab never actually
        renders survey geometry for MARE2DEM results loaded this way.
        """
        with pytest.raises(AttributeError, match="utm"):
            ivr._make_fig_mare2dem(mare_result, "survey", {})

    def test_survey_tab_renders_when_survey_layout_succeeds(
        self, mare_result, monkeypatch
    ):
        """Complements the "broken" test above: with a working
        PlotSurveyLayout, the callback's own fig/tight_layout/return
        plumbing for the survey tab is exercised."""
        import pycsamt.models.mare2dem.plot as mare_plot

        class _FakeSurveyLayout:
            def __init__(self, *a, **k):
                pass

            def plot(self, ax=None):
                pass

        monkeypatch.setattr(mare_plot, "PlotSurveyLayout", _FakeSurveyLayout)
        fig = ivr._make_fig_mare2dem(mare_result, "survey", {})
        assert fig is not None

    def test_unknown_tab_raises(self, mare_result):
        with pytest.raises(ValueError, match="Unknown MARE2DEM tab"):
            ivr._make_fig_mare2dem(mare_result, "not-a-tab", {})


# ── Callbacks: folder browser ─────────────────────────────────────────────


class TestBrowserToggle:
    def _fn(self, web_app):
        return _cb_multi(web_app, "invr-browser-modal.is_open")

    def test_browse_opens_at_home_by_default(self, web_app, monkeypatch):
        import pycsamt.app.web.callbacks.inv_results as ivr_mod

        monkeypatch.setattr(
            ivr_mod,
            "ctx",
            type("C", (), {"triggered_id": "btn-invr-browse"})(),
        )
        is_open, path = self._fn(web_app)(1, 0, False, None, None)
        assert is_open is True

    def test_browse_uses_typed_path_dir(self, web_app, monkeypatch, tmp_path):
        import pycsamt.app.web.callbacks.inv_results as ivr_mod

        monkeypatch.setattr(
            ivr_mod,
            "ctx",
            type("C", (), {"triggered_id": "btn-invr-browse"})(),
        )
        is_open, path = self._fn(web_app)(1, 0, False, None, str(tmp_path))
        assert is_open is True
        assert path == str(tmp_path)

    def test_cancel_closes(self, web_app, monkeypatch):
        import pycsamt.app.web.callbacks.inv_results as ivr_mod

        monkeypatch.setattr(
            ivr_mod,
            "ctx",
            type("C", (), {"triggered_id": "btn-invr-browser-cancel"})(),
        )
        is_open, path = self._fn(web_app)(0, 1, True, "/some/path", None)
        assert is_open is False
        assert path == "/some/path"

    def test_browse_falls_back_to_cur_path_when_no_typed_path(
        self, web_app, monkeypatch, tmp_path
    ):
        import pycsamt.app.web.callbacks.inv_results as ivr_mod

        monkeypatch.setattr(
            ivr_mod,
            "ctx",
            type("C", (), {"triggered_id": "btn-invr-browse"})(),
        )
        is_open, path = self._fn(web_app)(1, 0, False, str(tmp_path), None)
        assert is_open is True
        assert path == str(tmp_path)

    def test_browse_with_invalid_typed_path_falls_back_to_home(
        self, web_app, monkeypatch
    ):
        import os

        import pycsamt.app.web.callbacks.inv_results as ivr_mod

        monkeypatch.setattr(
            ivr_mod,
            "ctx",
            type("C", (), {"triggered_id": "btn-invr-browse"})(),
        )
        is_open, path = self._fn(web_app)(
            1, 0, False, None, "/no/such/typed/path/xyz"
        )
        assert is_open is True
        assert path == os.path.expanduser("~")

    def test_browse_typed_path_is_a_file_uses_its_parent_dir(
        self, web_app, monkeypatch, tmp_path
    ):
        import pycsamt.app.web.callbacks.inv_results as ivr_mod

        f = tmp_path / "somefile.txt"
        f.write_text("x")
        monkeypatch.setattr(
            ivr_mod,
            "ctx",
            type("C", (), {"triggered_id": "btn-invr-browse"})(),
        )
        is_open, path = self._fn(web_app)(1, 0, False, None, str(f))
        assert is_open is True
        assert path == str(tmp_path)

    def test_browse_with_zero_n_clicks_returns_closed(self, web_app, monkeypatch):
        import pycsamt.app.web.callbacks.inv_results as ivr_mod

        monkeypatch.setattr(
            ivr_mod,
            "ctx",
            type("C", (), {"triggered_id": "btn-invr-browse"})(),
        )
        is_open, path = self._fn(web_app)(0, 0, True, "/some/path", None)
        assert is_open is False
        assert path == "/some/path"


class TestBrowserNavigate:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "invr-browser-path.data", "invr-dir-item")

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)([None, None])

    def test_triggered_id_not_a_dict_prevents_update(self, web_app, monkeypatch):
        import pycsamt.app.web.callbacks.inv_results as ivr_mod

        monkeypatch.setattr(
            ivr_mod, "ctx", type("C", (), {"triggered_id": None})()
        )
        with pytest.raises(PreventUpdate):
            self._fn(web_app)([1])

    def test_navigates_to_clicked_dir(self, web_app, monkeypatch, tmp_path):
        import pycsamt.app.web.callbacks.inv_results as ivr_mod

        sub = tmp_path / "sub"
        sub.mkdir()
        monkeypatch.setattr(
            ivr_mod,
            "ctx",
            type(
                "C",
                (),
                {"triggered_id": {"type": "invr-dir-item", "index": str(sub)}},
            )(),
        )
        out = self._fn(web_app)([1])
        assert out == str(sub)

    def test_nonexistent_target_prevents_update(self, web_app, monkeypatch):
        import pycsamt.app.web.callbacks.inv_results as ivr_mod

        monkeypatch.setattr(
            ivr_mod,
            "ctx",
            type(
                "C",
                (),
                {
                    "triggered_id": {
                        "type": "invr-dir-item",
                        "index": "/no/such/dir",
                    }
                },
            )(),
        )
        with pytest.raises(PreventUpdate):
            self._fn(web_app)([1])


class TestBrowserUp:
    def _fn(self, web_app):
        return _cb_by_input(web_app, "invr-browser-path.data", "btn-invr-browser-up")

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "/a/b")

    def test_no_cur_path_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, None)

    def test_goes_up_one_level(self, web_app):
        out = self._fn(web_app)(1, str(Path("a") / "b"))
        assert out == "a"

    def test_at_root_prevents_update(self, web_app):
        root = Path(Path.cwd().anchor)
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, str(root))


class TestBrowserRefresh:
    def _fn(self, web_app):
        return _cb_multi(web_app, "invr-browser-list.children")

    def test_uses_home_when_no_path(self, web_app):
        import os

        items, cwd = self._fn(web_app)(None)
        assert cwd == os.path.expanduser("~")

    def test_lists_given_path(self, web_app, tmp_path):
        (tmp_path / "child").mkdir()
        items, cwd = self._fn(web_app)(str(tmp_path))
        assert cwd == str(tmp_path)


class TestBrowserSelect:
    def _fn(self, web_app):
        return _cb_by_input(
            web_app, f"{IDs.INVR_PATH_INPUT}.value", "btn-invr-browser-select"
        )

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "/a/b")

    def test_selects_path_and_closes_modal(self, web_app):
        value, is_open = self._fn(web_app)(1, "/a/b")
        assert value == "/a/b"
        assert is_open is False


# ── Callback: load ────────────────────────────────────────────────────────


class TestLoadResult:
    def _fn(self, web_app):
        return _cb_multi(web_app, f"{IDs.INVR_STORE}.data")

    def test_no_clicks_or_no_path_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(None, "/a/b", "auto", "s")
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, "", "auto", "s")

    def test_path_not_found(self, web_app):
        out = self._fn(web_app)(1, "/no/such/dir/xyz", "auto", "s")
        assert "Path not found" in out[1][1]

    def test_unknown_solver(self, web_app, tmp_path):
        out = self._fn(web_app)(1, str(tmp_path), "auto", "s")
        assert "Cannot detect solver" in out[1][1]

    @pytest.mark.skipif(not _HAS_OCCAM, reason="occam2D data not available")
    def test_loads_occam2d_successfully(self, web_app):
        meta, status, info, opts, tab, spinner_style = self._fn(web_app)(
            1, str(_OCCAM_DIR), "auto", "invr-sess-occam"
        )
        assert meta["solver"] == "occam2d"
        assert tab == "conv"
        assert spinner_style == {"display": "none"}
        assert ivr._fetch_result("invr-sess-occam") is not None

    @pytest.mark.skipif(not _HAS_MARE, reason="MARE2DEM data not available")
    def test_loads_mare2dem_successfully(self, web_app):
        meta, status, info, opts, tab, spinner_style = self._fn(web_app)(
            1, str(_MARE_DIR), "mare2dem", "invr-sess-mare"
        )
        assert meta["solver"] == "mare2dem"

    @pytest.mark.skipif(not _HAS_OCCAM, reason="occam2D data not available")
    def test_load_error_reported(self, web_app, monkeypatch):
        def _boom(path):
            raise RuntimeError("load boom")

        monkeypatch.setattr(ivr, "_LOADERS", {**ivr._LOADERS, "occam2d": _boom})
        out = self._fn(web_app)(1, str(_OCCAM_DIR), "occam2d", "s")
        assert "Load error" in out[1][1]

    def test_loads_modem_successfully_via_fake_loader(
        self, web_app, monkeypatch, tmp_path
    ):
        """ModEM sample data isn't shipped; stub the loader itself to
        exercise _load_result's modem success branch (station-options
        building, info strip, meta/tab construction)."""
        fake_result = _FakeModemResult()
        monkeypatch.setattr(
            ivr, "_LOADERS", {**ivr._LOADERS, "modem": lambda path: fake_result}
        )
        meta, status, info, opts, tab, spinner_style = self._fn(web_app)(
            1, str(tmp_path), "modem", "invr-sess-modem"
        )
        assert meta["solver"] == "modem"
        assert meta["n_iter"] == 10
        assert opts == [
            {"label": "mt01", "value": "mt01"},
            {"label": "mt02", "value": "mt02"},
        ]
        assert tab == "conv"
        assert spinner_style == {"display": "none"}
        texts = [str(getattr(i, "children", i)) for i in info]
        assert any("MODEM" in t for t in texts)
        assert ivr._fetch_result("invr-sess-modem") is fake_result

    def test_rms_conversion_error_falls_back_to_nan_and_question_mark(
        self, web_app, monkeypatch, tmp_path
    ):
        fake_result = _FakeModemResult()
        fake_result.final_rms = "not-a-number"
        monkeypatch.setattr(
            ivr, "_LOADERS", {**ivr._LOADERS, "modem": lambda path: fake_result}
        )
        meta, status, info, opts, tab, spinner_style = self._fn(web_app)(
            1, str(tmp_path), "modem", "invr-sess-badrms"
        )
        assert meta["rms"] != meta["rms"]  # NaN
        assert "RMS=?" in status[1]

    def test_loads_modem_without_data_station_opts_stay_empty(
        self, web_app, monkeypatch, tmp_path
    ):
        fake_result = _FakeModemResult()
        fake_result.data_obs = None
        fake_result.data_pred = None
        monkeypatch.setattr(
            ivr, "_LOADERS", {**ivr._LOADERS, "modem": lambda path: fake_result}
        )
        meta, status, info, opts, tab, spinner_style = self._fn(web_app)(
            1, str(tmp_path), "modem", "invr-sess-nodata"
        )
        assert meta["solver"] == "modem"
        assert opts == []

    def test_modem_station_opts_exception_is_swallowed(
        self, web_app, monkeypatch, tmp_path
    ):
        class _BadStations:
            @property
            def stations(self):
                raise RuntimeError("boom")

        fake_result = _FakeModemResult()
        fake_result.data_obs = _BadStations()
        fake_result.data_pred = None
        monkeypatch.setattr(
            ivr, "_LOADERS", {**ivr._LOADERS, "modem": lambda path: fake_result}
        )
        meta, status, info, opts, tab, spinner_style = self._fn(web_app)(
            1, str(tmp_path), "modem", "invr-sess-badstations"
        )
        assert meta["solver"] == "modem"  # load still succeeds
        assert opts == []

    def test_loads_occam2d_without_response_station_opts_stay_empty(
        self, web_app, monkeypatch, tmp_path
    ):
        class _FakeOccamResult:
            response = None
            n_iterations = 3
            final_rms = 2.0

        monkeypatch.setattr(
            ivr,
            "_LOADERS",
            {**ivr._LOADERS, "occam2d": lambda path: _FakeOccamResult()},
        )
        meta, status, info, opts, tab, spinner_style = self._fn(web_app)(
            1, str(tmp_path), "occam2d", "invr-sess-occam-noresp"
        )
        assert meta["solver"] == "occam2d"
        assert opts == []


# ── Callback: generate plot ────────────────────────────────────────────────


class TestGeneratePlot:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.INVR_CANVAS}.src")

    def test_no_meta_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, "conv", None, "s", *_default_ctrl_args())

    def test_no_cached_result_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                1,
                "conv",
                {"solver": "occam2d"},
                "no-such-session",
                *_default_ctrl_args(),
            )

    def test_unknown_solver_prevents_update(self, web_app):
        ivr._store_result("s1", object())
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                1,
                "conv",
                {"solver": "no-such-solver"},
                "s1",
                *_default_ctrl_args(),
            )

    @pytest.mark.skipif(not _HAS_OCCAM, reason="occam2D data not available")
    def test_real_occam2d_render(self, web_app, occam_result):
        ivr._store_result("s-occam", occam_result)
        src = self._fn(web_app)(
            1,
            "model",
            {"solver": "occam2d"},
            "s-occam",
            *_default_ctrl_args(),
        )
        assert src.startswith("data:image/png")

    @pytest.mark.skipif(not _HAS_OCCAM, reason="occam2D data not available")
    def test_plot_exception_falls_back_to_error_fig(self, web_app, occam_result):
        ivr._store_result("s-occam-err", occam_result)
        src = self._fn(web_app)(
            1,
            "not-a-real-tab",
            {"solver": "occam2d"},
            "s-occam-err",
            *_default_ctrl_args(),
        )
        assert src.startswith("data:image/png")

    def test_real_modem_render_via_fakes(
        self, web_app, fake_modem_result, fake_modem_plots
    ):
        ivr._store_result("s-modem", fake_modem_result)
        src = self._fn(web_app)(
            1,
            "conv",
            {"solver": "modem"},
            "s-modem",
            *_default_ctrl_args(),
        )
        assert src.startswith("data:image/png")


# ── Callback: export plot ──────────────────────────────────────────────────


class TestExportPlot:
    def _fn(self, web_app):
        return _cb(web_app, f"{IDs.INVR_DOWNLOAD}.data")

    def test_no_clicks_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                None,
                {"solver": "occam2d"},
                "conv",
                "s",
                *_default_ctrl_args(),
            )

    def test_no_meta_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(1, None, "conv", "s", *_default_ctrl_args())

    def test_no_cached_result_prevents_update(self, web_app):
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                1,
                {"solver": "occam2d"},
                "conv",
                "no-such-session",
                *_default_ctrl_args(),
            )

    def test_unknown_solver_prevents_update(self, web_app):
        ivr._store_result("s2", object())
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                1,
                {"solver": "no-such-solver"},
                "conv",
                "s2",
                *_default_ctrl_args(),
            )

    @pytest.mark.skipif(not _HAS_OCCAM, reason="occam2D data not available")
    def test_real_export(self, web_app, occam_result):
        ivr._store_result("s-exp", occam_result)
        out = self._fn(web_app)(
            1,
            {"solver": "occam2d"},
            "model",
            "s-exp",
            *_default_ctrl_args(),
        )
        assert out["filename"] == "invr_occam2d_model.png"
        assert out["base64"] is True

    @pytest.mark.skipif(not _HAS_OCCAM, reason="occam2D data not available")
    def test_export_exception_prevents_update(self, web_app, occam_result):
        ivr._store_result("s-exp-err", occam_result)
        with pytest.raises(PreventUpdate):
            self._fn(web_app)(
                1,
                {"solver": "occam2d"},
                "not-a-real-tab",
                "s-exp-err",
                *_default_ctrl_args(),
            )

    def test_real_modem_export_via_fakes(
        self, web_app, fake_modem_result, fake_modem_plots
    ):
        ivr._store_result("s-modem-exp", fake_modem_result)
        out = self._fn(web_app)(
            1,
            {"solver": "modem"},
            "conv",
            "s-modem-exp",
            *_default_ctrl_args(),
        )
        assert out["filename"] == "invr_modem_conv.png"
        assert out["base64"] is True
